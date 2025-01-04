
export demo
export spanwise_sln

function xcenter(quad)
    return 0.25 * sum(map(vert -> coords(vert).x, vertices(quad)))
end

# Only camber, not thickness included!
function nacacamberlaw(nacanums, pos)
    # Simply derivative of thickness function
    t = (nacanums[3] * 10 + nacanums[4]) / 100
    thickprime = 5t * (-0.126 + 0.14845 / sqrt(pos) - 0.7032pos + 0.8529pos^2 - 0.406pos^3)

    # note, derivative of the position is dy/dx, thus angle of slope is 
    # the arctangent of this value!
    if nacanums[1] == 0 && nacanums[2] == 0
        return 0.0
    else
        m = nacanums[1] / 100
        p = nacanums[2] / 10
        camberprime = 0.0
        if pos <= p
            camberprime = m / p^2 * (2p - 2pos)
        else
            camberprime = m / (1 - p)^2 * (2p - 2pos)
        end

        return atan(camberprime)
    end
end

# Ad-hoc and inefficient. 
function spanwise_sln(mesh, sln)
    X = []
    Y = []

    for i in eachindex(mesh)
        found = false
        nx = xcenter(mesh[i])
        for xi in eachindex(X)
            if X[xi] == nx
                Y[xi] += sln[i]
                found = true
                break
            end
        end
        if !found
            push!(X, nx)
            push!(Y, sln[i])
        end
    end

    return ustrip(hcat(X, Y))
end

function visualizescalar(comb, slopes)
    display(viz(comb, color=slopes))
end

function demo(α)
    grid1, slopes1 = wingmesh(40, 40, 15, 1.9, 0.8, deg2rad(15), x -> nacacamberlaw([4, 4, 1, 6], x))
    grid2, slopes2 = wingmesh(30, 30, 7.5, 0.9, 0.6, deg2rad(10), x -> nacacamberlaw([0, 0, 1, 8], x))
    grid2 = Translate(0, 7)(grid2)
    comb = grid1 ∪ grid2

    slopes = vcat(slopes1, slopes2)

    θs = zeros(length(comb))

    # Only first wing imposes theta on its solution
    for (i, elem) in enumerate(grid1)
        x = xcenter(elem)
        xadim = x / 15
        θs[i] = 2.0 * abs(ustrip(xadim)) * deg2rad(5)
    end

    rh = rhs(mesh, α, slopes, θs)
    #visualizescalar(comb, rh)

    sln = solve(comb, α, slopes, θs)

    #display(viz(comb, color=log.(abs.(sln))))
    return coefficients(comb, sln, α, 0)
end
