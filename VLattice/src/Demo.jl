
export demo
export spanwise_sln

function xcenter(quad)
    return 0.25 * sum(map(vert -> coords(vert).x, vertices(quad)))
end

# Ad-hoc and inefficient
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

function demo(α)
    grid1 = wingmesh(50, 50, 15, 1.9, 0.8, deg2rad(15))
    grid2 = wingmesh(40, 40, 7.5, 0.9, 0.6, deg2rad(10))
    grid2 = Translate(0, 7)(grid2)
    comb = grid1 ∪ grid2

    θs = zeros(length(comb))

    # Only first wing imposes theta on its solution
    for (i, elem) in enumerate(grid1)
        x = xcenter(elem)
        xadim = x / 15
        θs[i] = abs(ustrip(xadim)) * deg2rad(5)
    end

    sln = solve(comb, α, θs)
    return comb, sln
    #display(viz(comb, color=log.(abs.(sln))))
end
