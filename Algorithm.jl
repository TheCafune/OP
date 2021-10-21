include("Graph.jl")
include("Matrix.jl")

using LinearAlgebra
using SparseArrays
using Laplacians

function getResult(G, ff, k, S, s)
    sf = getSf(G.n, k, s, ff)
    ss = getSs(G.n, k, s, ff)
    Omg = getOmegaS(G, k, ff, S)
    Afs = getAfsS(G, k, ff, S)
    of = ones(G.n-k)
    ans = of' * Omg * (Afs * ss + sf)
    return ans
end

function randomSelect(Ev, k2; randomSeed = round(Int, time() * 10000))
    Random.seed!(randomSeed)
    y = randperm(size(Ev, 1))
    S = Array{Tuple{Int32, Int32}, 1}()
    for i = 1 : k2
        push!(S, Ev[y[i]])
    end
    return S
end

function exactOpinion(G, ff, k, k2, s, Ev)
    sf = getSf(G.n, k, s, ff)
    ss = getSs(G.n, k, s, ff)
    Omg = getOmega(G, k, ff)
    Afs = getAfs(G, k, ff)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    of = ones(G.n-k)
    for i = 1 : k2
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = (1.0 / (1.0+Omg[uu, uu])) * (of' * Omg[:, uu]) * (1.0 - Omg[uu, :]' * (Afs * ss + sf))
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        (su, sv) = Ev[xx]
        Omg = Omg - (1.0 / (1.0 + Omg[su, su])) * Omg[:, su] * Omg[su, :]'
        Afs[su, sv-G.n+k] += 1.0
        push!(S, Ev[xx])
    end
    return S
end

function approxOpinion(G, ff, k, k2, s, Ev; eps = 0.3)
    sf = getSf(G.n, k, s, ff)
    ss = getSs(G.n, k, s, ff)
    IpL = getSparseIpL(G, k, ff)
    sAfs = getSparseAfs(G, k, ff)
    B, X = getSparseBandX(G, k, ff)
    m0 = size(B, 1)
    kkk = round(Int, log2(G.n) / (eps^2))
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    of = ones(G.n-k)
    for rep = 1 : k2
        diag = zeros(G.n-k)
        f = approxchol_sddm(IpL, tol=1e-12)
        for i = 1 : kkk
            yy1 = B' * randn(m0)
            yy2 = X * randn(G.n-k)
            zz1 = f(yy1)
            zz2 = f(yy2)
            for j = 1 : G.n-k
                diag[j] += (zz1[j]^2 + zz2[j]^2)
            end
        end
        diag ./= kkk
        of2 = f(of)
        of3 = f(sAfs * ss + sf)

        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = (1.0 / (1.0+diag[uu])) * (of2[uu]) * (1.0 - of3[uu])
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        (su, sv) = Ev[xx]
        IpL[su, su] += 1
        sAfs[su, sv-G.n+k] += 1.0
        push!(S, Ev[xx])
    end
    return S
end

function topDegree(G, ff, k, k2, Ev)
    d = zeros(G.n)
    for (ID, u, v, w) in G.E
        nu = ff[u]
        nv = ff[v]
        d[nu] += w
        d[nv] += w
    end
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for rep = 1 : k2
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                dz[j] = -10000
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = d[uu]
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        (su, sv) = Ev[xx]
        d[su] += 1
        push!(S, Ev[xx])
    end
    return S
end

function BetweennessCentrality(G, ff, S)
    gg = zeros(Int, G.n)
    foreach(i -> gg[ff[i]] = i, 1 : G.n)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    for (u, v) in S
        push!(g[gg[u]], gg[v])
        push!(g[gg[v]], gg[u])
    end
    C = zeros(G.n)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    Q = zeros(Int32, G.n+10)
    delta = zeros(G.n)
    for s = 1 : G.n
        foreach(i -> p[i] = [], 1 : G.n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                if w != s
                    C[w] += delta[w]
                end
            end
        end

    end

    return C
end

function topBetweenness(G, ff, k, k2, Ev)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for rep = 1 : k2
        dd = BetweennessCentrality(G, ff, S)
        d = zeros(G.n)
        foreach(i -> d[ff[i]] = dd[i], 1 : G.n)
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                dz[j] = -10000
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = d[uu]
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        push!(S, Ev[xx])
    end
    return S
end

function alpP(G, ff, S, alp)
    gg = zeros(Int, G.n)
    foreach(i -> gg[ff[i]] = i, 1 : G.n)
    d = zeros(G.n)
    for (ID, u, v, w) in G.E
        d[u] += w
        d[v] += w
    end
    for (u, v) in S
        d[gg[u]] += 1
        d[gg[v]] += 1
    end
    mk = size(S, 1)
    Is = zeros(Int, G.m*2 + mk*2)
    Js = zeros(Int, G.m*2 + mk*2)
    Vs = zeros(G.m*2 + mk*2)
    for (ID, u, v, w) in G.E
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = alp * w / d[v]
        Is[G.m + ID] = v
        Js[G.m + ID] = u
        Vs[G.m + ID] = alp * w / d[u]
    end
    for i = 1 : mk
        (uu, vv) = S[i]
        u = gg[uu]
        v = gg[vv]
        Is[G.m*2 + i] = u
        Js[G.m*2 + i] = v
        Vs[G.m*2 + i] = alp * 1.0 / d[v]
        Is[G.m*2 + mk + i] = v
        Js[G.m*2 + mk + i] = u
        Vs[G.m*2 + mk + i] = alp * 1.0 / d[u]
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end

function PageRank(G, ff, S; alpha = 0.85)
    aP = alpP(G, ff, S, alpha)
    C = zeros(G.n)
    C[1] = 1.0
    adC = ((1.0 - alpha) / G.n) * ones(G.n)

    while true
        pC = copy(C)
        C = aP * C + adC
        if norm(C - pC) < 1e-12
            break
        end
    end

    return C
end

function topPageRank(G, ff, k, k2, Ev)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for rep = 1 : k2
        dd = PageRank(G, ff, S)
        d = zeros(G.n)
        foreach(i -> d[ff[i]] = dd[i], 1 : G.n)
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                dz[j] = -10000
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = d[uu]
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        push!(S, Ev[xx])
    end
    return S
end

function Opt(G, ff, k, k2, s, Ev)
    sf = getSf(G.n, k, s, ff)
    ss = getSs(G.n, k, s, ff)
    Omg = getOmega(G, k, ff)
    Afs = getAfs(G, k, ff)
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    of = ones(G.n-k)
    bst = 0.0

    optsn = zeros(Int, k2)
    a = zeros(Int, k2)
    stkOmg = []
    stkAfs = []
    push!(stkOmg, Omg)
    push!(stkAfs, Afs)
    top = 1

    dfs(dep) = begin
        if dep > k2
            tmp = of' * stkOmg[top] * (stkAfs[top] * ss + sf)
            if tmp > bst
                bst = tmp
                for i = 1 : k2
                    optsn[i] = a[i]
                end
            end
            return nothing
        end
        st = 1
        if dep > 1
            st = a[dep-1]
        end
        for i = st : sev
            if !cho[i]
                cho[i] = true
                (su, sv) = Ev[i]
                nOmg = copy(stkOmg[top])
                nOmg = nOmg - (1.0 / (1.0 + nOmg[su, su])) * nOmg[:, su] * nOmg[su, :]'
                nAfs = copy(stkAfs[top])
                nAfs[su, sv-G.n+k] += 1.0
                top += 1
                push!(stkOmg, nOmg)
                push!(stkAfs, nAfs)
                a[dep] = i
                dfs(dep+1)
                pop!(stkOmg)
                pop!(stkAfs)
                top -= 1
                cho[i] = false
            end
        end
        return nothing
    end

    dfs(1)
    S = Array{Tuple{Int32, Int32}, 1}()
    for i = 1 : k2
        push!(S, Ev[optsn[i]])
    end
    return S
end

function ClosenessCentrality(G, ff, S)
    gg = zeros(Int, G.n)
    foreach(i -> gg[ff[i]] = i, 1 : G.n)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    for (u, v) in S
        push!(g[gg[u]], gg[v])
        push!(g[gg[v]], gg[u])
    end
    C = zeros(G.n)
    d = zeros(Int32, G.n)
    Q = zeros(Int32, G.n+10)
    for s = 1 : G.n
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
            end
        end

        C[s] = sum(d)
    end

    foreach(i -> C[i] = 1.0 / C[i], 1 : G.n)

    return C
end

function topCloseness(G, ff, k, k2, Ev)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for rep = 1 : k2
        dd = ClosenessCentrality(G, ff, S)
        d = zeros(G.n)
        foreach(i -> d[ff[i]] = dd[i], 1 : G.n)
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                dz[j] = -10000
                continue
            end
            (uu, vv) = Ev[j]
            dz[j] = d[uu]
        end
        xx = argmax(dz)
        #println(dz[xx])
        cho[xx] = true
        push!(S, Ev[xx])
    end
    return S
end
