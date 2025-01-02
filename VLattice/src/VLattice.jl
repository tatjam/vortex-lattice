module VLattice

using Meshes
using Unitful
using GLMakie

export wingmesh
export quad2elem
export influence
export influencemat
export solve

"""
	wingmesh(nspanwise, nchordwise, b, Γ, claw, cofflaw)

Generates a structured mesh of given vertex numbers for a wing of span `b`, and with the given 
functions `claw` and `cofflaw` which define the distribution of chord and local 1/4 point offset at 
each spanwise position.
These are functions of dimensional span position, going from -b/2 to b/2.

The wing span-wise direction is assumed to be parallel to the x-axis.
If cofflaw = 0, the 1/4 point of the chord will lie on the X axis.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, claw::Function, cofflaw::Function)
    X = zeros(nchordwise, nspanwise)
    Y = zeros(nchordwise, nspanwise)
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
        end
    end

    return StructuredGrid(X, Y)
end

"""
	wingmesh(nspawnwise, nchordwise, b, claw, Γ)

Same as generic `wingmesh` but generates the wing given its sweep in radians. Sweep is the angle 
formed by the 1/4 chord points of the root and the tip of the wing.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, claw::Function, Γ::Number)
    return wingmesh(nspanwise, nchordwise, b, claw, x -> abs(x) * tan(Γ))
end

"""
	wingmesh(nspawnwise, nchordwise, rootchord, tipchord, Γ)

Same as generic `wingmesh` but generates the wing given its sweep in radians. Sweep is the angle 
formed by the 1/4 chord points of the root and the tip of the wing.
"""
function wingmesh(nspanwise::Integer, nchordwise::Integer, b::Number, rootchord::Number,
    tipchord::Number, Γ::Number)
    return wingmesh(nspanwise, nchordwise, b,
        x -> rootchord - (rootchord - tipchord) * abs(x) / (0.5 * b), Γ)
end

function leading_trailing(pa, pb)
    if coords(pa).y > coords(pb).y
        return [pa, pb]
    else
        return [pb, pa]
    end
end

function testpoint(quad::Quadrangle)
    elem = vertices(quad)

    # 1/4 point in y direction of central x point. This doesn't matter that much actually.	
    xmiddle = 0.5 * (coords(elem[1]).x + coords(elem[4]).x)
    ylead = 0.5 * (coords(elem[1]).y + coords(elem[4]).y)
    ytrail = 0.5 * (coords(elem[2]).y + coords(elem[3]).y)

    return [ustrip(xmiddle), ustrip(ylead * 0.75 + ytrail * 0.25)]
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
    celem = vertices(cause)

    # From class notes, but changing x for y
    a = tp[2] - ustrip(coords(celem[1]).y)
    b = tp[1] - ustrip(coords(celem[1]).x)
    c = tp[2] - ustrip(coords(celem[4]).y)
    d = tp[1] - ustrip(coords(celem[4]).x)

    e = sqrt(a^2 + b^2)
    f = sqrt(c^2 + d^2)

    g = ustrip(coords(celem[4]).y) - ustrip(coords(celem[1]).y)
    h = ustrip(coords(celem[4]).x) - ustrip(coords(celem[1]).x)

    val1 = (g * a + h * b) / e - (g * c + h * d) / f
    val2 = (c + f) / (d * f) - (a + e) / (b * e)
    val = val1 / (a * d - c * b) + val2

    return val / (4 * π)
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

We impose the condition that the flow is parallel to the tilted plane, which 
roughly means that, taking `Vlocal = V∞ + Vinduced`, `tan(Vlocal_z / Vlocal_y)`
is roughly equivalent to `Vinduced / V∞`, and this value must be equal to alpha.
We set V∞ to 1, so the equation turns out to be Vinduced = alpha, so the right 
hand side is simply alpha

"""
function rhs(mesh, α)
    totalsize = length(mesh)
    return fill(α, totalsize)
end

function solve(mesh, α)
    mat = influencemat(mesh)
    rh = rhs(mesh, α)
    return mat \ rh
end

end # module VLattice
