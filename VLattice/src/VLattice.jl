module VLattice

using Meshes
using GLMakie

export wingmesh

"""
	wingmesh(nspanwise, nchordwise, b, Γ, claw, cofflaw)

Generates a structured mesh of given vertex numbers for a wing of span `b`, and with the given 
functions `claw` and `cofflaw` which define the distribution of chord and local 1/4 point offset at 
each spanwise position.
These are functions of dimensional span position, going from -b/2 to b/2.

The wing span-wise direction is assumed to be parallel to the x-axis.
If cofflaw = 0, the 1/4 point of the chord will lie on the X axis.
"""
function wingmesh(nspanwise, nchordwise, b, claw::Function, cofflaw::Function)
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
            chordpos = c * (yprog - 0.5) - qpos + coff
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
function wingmesh(nspanwise, nchordwise, b, claw::Function, Γ)
    return wingmesh(nspanwise, nchordwise, b, claw, x -> -abs(x) * tan(Γ))
end

end # module VLattice
