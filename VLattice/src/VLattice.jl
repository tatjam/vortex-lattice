module VLattice

using Meshes
using Unitful
using GLMakie
using LinearAlgebra

export wingmesh
export quad2elem
export influence
export influencemat
export solve
export coefficients

include("Demo.jl")

"""
	wingmesh(nspanwise, nchordwise, b, Γ, claw, cofflaw, camberlaw)

Generates a structured mesh of given vertex numbers for a wing of span `b`, and with the given 
functions `claw` and `cofflaw` which define the distribution of chord and local 1/4 point offset at 
each spanwise position.
These are functions of dimensional span position, going from -b/2 to b/2. The camberlaw is 
instead given NON-DIMENSIONAL chord position, from 0 to 1

The wing span-wise direction is assumed to be parallel to the x-axis.
If cofflaw = 0, the 1/4 point of the chord will lie on the X axis.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, claw::Function,
    cofflaw::Function, camberlaw::Function)

    X = zeros(nchordwise, nspanwise)
    Y = zeros(nchordwise, nspanwise)
    slopes = zeros((nchordwise - 1) * (nspanwise - 1))

    for i in 1:nspanwise
        linprog = (i - 1) / (nspanwise - 1) # from 0 to 1
        # Adjust progression to accumulate vertices near the edges
        prog = 0.5 * (1 - cos(π * linprog)) # from 0 to 1

        spanpos = b * (prog - 0.5) # in units of length

        # Spanwise vertices are all parallel
        X[:, i] .= spanpos

        c = claw(spanpos)
        coff = cofflaw(spanpos)

        for yi in 1:nchordwise
            linyprog = (yi - 1) / (nchordwise - 1) # from 0 to 1
            yprog = 0.5 * (1 - cos(π * linyprog)) # from 0 to 1

            qpos = c * 0.25
            chordpos = c * (yprog - 0.5) + qpos + coff
            Y[yi, i] = chordpos

            if i > 1 && i < nspanwise
                if yi > 1 && yi < nchordwise
                    slopes[(i-1)*(nchordwise-1)+(yi-1)] = camberlaw(yprog)
                end
            end
        end
    end

    return StructuredGrid(X, Y), slopes
end

"""
	wingmesh(nspawnwise, nchordwise, b, claw, Γ)

Same as generic `wingmesh` but generates the wing given its sweep in radians. Sweep is the angle 
formed by the 1/4 chord points of the root and the tip of the wing.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, claw::Function, Γ::Number, camberlaw::Function)
    return wingmesh(nspanwise, nchordwise, b, claw, x -> abs(x) * tan(Γ), camberlaw)
end

"""
	wingmesh(nspawnwise, nchordwise, rootchord, tipchord, Γ, camberlaw)

Same as generic `wingmesh` but generates the wing given its sweep in radians. Sweep is the angle 
formed by the 1/4 chord points of the root and the tip of the wing.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, rootchord::Number,
    tipchord::Number, Γ::Number, camberlaw::Function)
    return wingmesh(nspanwise, nchordwise, b,
        x -> rootchord - (rootchord - tipchord) * abs(x) / (0.5 * b), Γ, camberlaw)
end

function testpoint(quad::Quadrangle)
    elem = vertices(quad)

    # 1/4 point in y direction of central x point. This doesn't matter that much actually.	
    xmiddle = 0.5 * (coords(elem[1]).x + coords(elem[4]).x)
    ylead = 0.5 * (coords(elem[1]).y + coords(elem[4]).y)
    ytrail = 0.5 * (coords(elem[2]).y + coords(elem[3]).y)

    return [ustrip(xmiddle), ustrip(ylead * 0.25 + ytrail * 0.75)]
end

function vortexpoints(quad::Quadrangle)
    elem = vertices(quad)

    return [
        [ustrip(coords(elem[1]).x), ustrip(coords(elem[1]).y * 0.75 + coords(elem[2]).y * 0.25)],
        [ustrip(coords(elem[4]).x), ustrip(coords(elem[4]).y * 0.75 + coords(elem[3]).y * 0.25)]
    ]
end

function vortexpoint(quad::Quadrangle)
    vpts = vortexpoints(quad)

    return [
        0.5 * (vpts[1][1] + vpts[2][1]),
        0.5 * (vpts[1][2] + vpts[2][2])
    ]
end

function levec(quad::Quadrangle)
    elem = vertices(quad)

    Δx = coords(elem[4]).x - coords(elem[1]).x
    Δy = coords(elem[4]).y - coords(elem[1]).y

    return [ustrip(Δx), ustrip(Δy), 0]
end

function lelen(quad::Quadrangle)
    elem = vertices(quad)

    Δx = coords(elem[4]).x - coords(elem[1]).x
    Δy = coords(elem[4]).y - coords(elem[1]).y
    l = sqrt(Δx * Δx + Δy * Δy)

    return ustrip(l)
end

"""
	influence(cause, tp)

Returns the induced velocity by a unitary horseshoe vertex shed by the given element of the grid,
on a given 2D position. Because everything is flat, this velocity is in the z-direction.
The vortex shedding direction is given by the y direction.
The leading edge is formed by the two elements with lesser y out of the two pairs with same x.

Assumes vertice are ordered as follows:
- Leading edge, trailing edge, trailing edge, leading edge
"""
function influence(cause::Quadrangle, tp)
    celem = vortexpoints(cause)

    # From class notes, but changing x for y
    a = tp[2] - celem[1][2]
    b = tp[1] - celem[1][1]
    c = tp[2] - celem[2][2]
    d = tp[1] - celem[2][1]

    e = sqrt(a^2 + b^2)
    f = sqrt(c^2 + d^2)

    g = celem[2][2] - celem[1][2]
    h = celem[2][1] - celem[1][1]

    val1 = (g * a + h * b) / e - (g * c + h * d) / f
    val2 = (c + f) / (d * f) - (a + e) / (b * e)
    val = val1 / (a * d - c * b) + val2

    if isnan(val) || isinf(val)
        return 0.0
    end

    return val / (4.0 * π)
end

"""
	influencemat(mesh)

Generates the influence matrix of every element of `mesh` on each other and with itself. The 
indexing of the resulting matrix is [effect, cause]

Use `union(mesh1, mesh2)` to join two meshes. If a single mesh is used, call with elements(mesh).

"""
function influencemat(mesh)
    totalsize = length(mesh)
    out = zeros(totalsize, totalsize)

    tp = zeros(totalsize, 2)
    # Generate test points for each element. This saves quite a bit of computation time
    for (i, panel) in enumerate(mesh)
        tp[i, 1:2] = testpoint(panel)
    end

    for (causei, cause) in enumerate(mesh)
        for effecti in eachindex(mesh)
            out[effecti, causei] += influence(cause, tp[effecti, :])
        end
    end

    return out
end

"""

We impose the condition that the flow is parallel to the tilted plane, i.e. 
vlocal⋅norm = 0 

where   vlocal = (0, cos(α), vind-sin(α))
        norm = (0, sin(θ + slope), cos(θ + slope))

thus vlocal⋅norm = sin(θ + slope)cos(α) + cos(θ + slope)*(vind - sin(α)) = 0
we can reorder terms:
    sin(α) - tan(θ + slope)cos(α) = vind
where vind is each of our matrix rows

"""
function rhs(mesh, α, slopes, θs=nothing)
    if isnothing(θs)
        θs = zeros(length(mesh))
    end
    return (sin(α) .- tan.(θs .- slopes) .* cos(α))
end

function solve(mesh, α, slopes, θs=nothing)
    mat = influencemat(mesh)
    rh = rhs(mesh, α, slopes, θs)
    return mat \ rh
end

# returns cL, cDi, pitch moment (with respect to given point)
# Assume 1m/s freestream velocity
function coefficients(mesh, sln, α, moment0)
    tot_surface = 0.0
    force = [0.0, 0.0, 0.0]
    forces = zeros(length(mesh))
    for (i, elem) in enumerate(mesh)
        s = ustrip(area(elem))
        tot_surface += s

        indvel = 0
        for (otheri, otherelem) in enumerate(mesh)
            if otheri != i
                indvel += influence(otherelem, vortexpoint(elem)) * sln[otheri]
            end
        end

        Vinf = [0, cos(α), indvel + sin(α)]
        # This is really force / rho 
        force .+= sln[i] .* cross(Vinf, levec(elem))
        forces[i] = sln[i] * cross(Vinf, levec(elem))[3]
    end

    # project force
    lift = force[2] * sin(α) + force[3] * cos(α)
    drag = force[2] * cos(α) - force[3] * sin(α)

    # cL / cD is force / 0.5 * rho * U^2 * S, thus we have to divide by surface and multiply by 2
    return 2.0 * lift / tot_surface, 2.0 * drag / tot_surface

end

end # module VLattice
