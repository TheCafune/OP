include("Graph.jl")

using Random

function generateLeader(n, k)
    Random.seed!(round(Int, time() * 10000))
    x = randperm(n)
    ld = zeros(Int, n)
    for i = 1 : k
        ld[x[i]] = 1
    end
    return ld
end

function generateOpinion(n, ld)
    Random.seed!(round(Int, time() * 10000))
    s = zeros(n)
    for i = 1 : n
        if ld[i] == 1
            s[i] = 1.0
        else
            s[i] = rand()
        end
    end
    return s
end

function generateEv(G, ff, k, evn)
    Random.seed!(round(Int, time() * 10000))
    nE = Set{Tuple{Int32, Int32}}()
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) && (nv<=(G.n-k))
            push!(nE, (nv, nu))
        end
        if (nu<=(G.n-k)) && (nv>(G.n-k))
            push!(nE, (nu, nv))
        end
    end
    sz = evn
    Ev = Array{Tuple{Int32, Int32}, 1}()
    while size(Ev, 1) < sz
        ru = rand(1:G.n-k)
        rv = rand(G.n-k+1:G.n)
        if !((ru, rv) in nE)
            push!(Ev, (ru, rv))
            push!(nE, (ru, rv))
        end
    end
    return Ev
end
