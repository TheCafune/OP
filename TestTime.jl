include("Graph.jl")
include("Algorithm.jl")
include("Gen.jl")
include("Matrix.jl")

function runMain(ags1, ags2)
    # Read Graph
    buf = split(ags1, ',')
    fileName = string("data/", buf[1], ".txt")
    networkType = "unweighted"
    if size(buf, 1) > 1
    	networkType = buf[2]
    end
    # Find LLC
    G0 = readGraph(fileName, networkType)
    G = getLLC(G0)

    buf2 = split(ags2, ',')
    k = parse(Int, buf2[1])
    mk2 = parse(Int, buf2[2])
    evn = parse(Int, buf2[3])

    # Generate Data
    ld = generateLeader(G.n, k)
    s = generateOpinion(G.n, ld)
    ff, gg = reLabel(G.n, ld)
    Ev = generateEv(G, ff, k, evn)

    lg = open("logTime.txt", "a")
    println(lg, buf[1])
    println(lg, G0.n, " ", G0.m, " ", G.n, " ", G.m)

    Time = zeros(2)
    if G.n < 50000
        Time[1] = time()
        S1 = exactOpinion(G, ff, k, mk2, s, Ev)
        Time[1] = time() - Time[1]
    end
    Time[2] = time()
    S2 = approxOpinion(G, ff, k, mk2, s, Ev)
    Time[2] = time() - Time[2]
    if G.n < 50000
        println(lg, "ExactOpinion : ", getResult(G, ff, k, S1, s))
        println(lg, "Time of ExactOpinion : ", Time[1])
    end
    println(lg, "ApproxOpinion : ", getResult(G, ff, k, S2, s))
    println(lg, "Time of ApproxOpinion : ", Time[2])
end

runMain(ARGS[1], ARGS[2])
