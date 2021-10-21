include("Graph.jl")

using SparseArrays

function reLabel(n, ld)
    ff = zeros(Int, n)
    gg = zeros(Int, n)
    cnt = 0
    cnt2 = n+1
    for i = 1 : n
        if ld[i] == 0
            cnt += 1
            ff[i] = cnt
            gg[cnt] = i
        else
            cnt2 -= 1
            ff[i] = cnt2
            gg[cnt2] = i
        end
    end
    return ff, gg
end

function getSf(n, k, s, ff)
    sf = zeros(n-k)
    for i = 1 : n
        if ff[i] <= (n-k)
            sf[ff[i]] = s[i]
        end
    end
    return sf
end

function getSs(n, k, s, ff)
    ss = zeros(k)
    for i = 1 : n
        if ff[i] > (n-k)
            ss[ff[i]-n+k] = s[i]
        end
    end
    return ss
end

function getOmega(G, k, ff)
    IpL = zeros(G.n-k, G.n-k)
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) || (nv>(G.n-k))
            if nu<=(G.n-k)
                IpL[nu, nu] += w
            end
            if nv<=(G.n-k)
                IpL[nv, nv] += w
            end
            continue
        end
        IpL[nu, nu] += w
        IpL[nv, nv] += w
        IpL[nu, nv] -= w
        IpL[nv, nu] -= w
    end
    for i = 1 : G.n-k
        IpL[i, i] += 1
    end
    return inv(IpL)
end

function getAfs(G, k, ff)
    Afs = zeros(G.n-k, k)
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) && (nv<=(G.n-k))
            Afs[nv, nu-G.n+k] = w
        end
        if (nu<=(G.n-k)) && (nv>(G.n-k))
            Afs[nu, nv-G.n+k] = w
        end
    end
    return Afs
end


function getOmegaS(G, k, ff, S)
    IpL = zeros(G.n-k, G.n-k)
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) || (nv>(G.n-k))
            if nu<=(G.n-k)
                IpL[nu, nu] += w
            end
            if nv<=(G.n-k)
                IpL[nv, nv] += w
            end
            continue
        end
        IpL[nu, nu] += w
        IpL[nv, nv] += w
        IpL[nu, nv] -= w
        IpL[nv, nu] -= w
    end
    for (u, v) in S
        IpL[u, u] += 1
    end
    for i = 1 : G.n-k
        IpL[i, i] += 1
    end
    return inv(IpL)
end

function getAfsS(G, k, ff, S)
    Afs = zeros(G.n-k, k)
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) && (nv<=(G.n-k))
            Afs[nv, nu-G.n+k] = w
        end
        if (nu<=(G.n-k)) && (nv>(G.n-k))
            Afs[nu, nv-G.n+k] = w
        end
    end
    for (u, v) in S
        Afs[u, v-G.n+k] = 1
    end
    return Afs
end

function getSparseIpL(G, k, ff)
    Is = Array{Int32, 1}()
    Js = Array{Int32, 1}()
    Vs = Array{Float64, 1}()
    for i = 1 : G.n-k
        push!(Is, i)
        push!(Js, i)
        push!(Vs, 1)
    end

    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) || (nv>(G.n-k))
            if nu<=(G.n-k)
                Vs[nu] += w
            end
            if nv<=(G.n-k)
                Vs[nv] += w
            end
            continue
        end
        Vs[nu] += w
        Vs[nv] += w
        push!(Is, nu)
        push!(Js, nv)
        push!(Vs, -w)
        push!(Is, nv)
        push!(Js, nu)
        push!(Vs, -w)
    end
    return sparse(Is, Js, Vs, G.n-k, G.n-k)
end

function getSparseAfs(G, k, ff)
    Is = Array{Int32, 1}()
    Js = Array{Int32, 1}()
    Vs = Array{Float64, 1}()
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) && (nv<=(G.n-k))
            push!(Is, nv)
            push!(Js, nu-G.n+k)
            push!(Vs, w)
        end
        if (nu<=(G.n-k)) && (nv>(G.n-k))
            push!(Is, nu)
            push!(Js, nv-G.n+k)
            push!(Vs, w)
        end
    end
    return sparse(Is, Js, Vs, G.n-k, k)
end

function getSparseBandX(G, k, ff)
    I1s = Array{Int32, 1}()
    J1s = Array{Int32, 1}()
    V1s = Array{Float64, 1}()
    I2s = zeros(Int32, G.n-k)
    J2s = zeros(Int32, G.n-k)
    V2s = ones(G.n-k)
    for i = 1 : G.n-k
        I2s[i] = i
        J2s[i] = i
    end
    m0 = 0
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        if (nu>(G.n-k)) || (nv>(G.n-k))
            if nu<=(G.n-k)
                V2s[nu] += w
            end
            if nv<=(G.n-k)
                V2s[nv] += w
            end
            continue
        end
        m0 += 1
        push!(I1s, m0)
        push!(J1s, nu)
        push!(V1s, sqrt(w))
        push!(I1s, m0)
        push!(J1s, nv)
        push!(V1s, -sqrt(w))
    end
    for i = 1 : G.n-k
        V2s[i] = sqrt(V2s[i])
    end
    B = sparse(I1s, J1s, V1s, m0, G.n-k)
    X = sparse(I2s, J2s, V2s, G.n-k, G.n-k)
    return B, X
end
