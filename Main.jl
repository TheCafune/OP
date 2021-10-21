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

    lg = open("log.txt", "a")
    println(lg, buf[1])
    println(lg, G0.n, " ", G0.m, " ", G.n, " ", G.m)

    rseed = round(Int, time() * 10000)
    for k2 = 1 : mk2
        println(lg, "k2 = ", k2)
        # Do Experiment
        S1 = randomSelect(Ev, k2; randomSeed = rseed)
        S2 = exactOpinion(G, ff, k, k2, s, Ev)
        S3 = approxOpinion(G, ff, k, k2, s, Ev)
        S4 = topDegree(G, ff, k, k2, Ev)
        if G.n < 600
            S5 = Opt(G, ff, k, k2, s, Ev)
        end
        S6 = topBetweenness(G, ff, k, k2, Ev)
        S7 = topPageRank(G, ff, k, k2, Ev)
        S8 = topCloseness(G, ff, k, k2, Ev)

        # Print Result
        println(lg, "Random : ", getResult(G, ff, k, S1, s))
        println(lg, "TopDegree : ", getResult(G, ff, k, S4, s))
        println(lg, "ExactOpinion : ", getResult(G, ff, k, S2, s))
        println(lg, "ApproxOpinion : ", getResult(G, ff, k, S3, s))
        if G.n < 600
            println(lg, "Optimum : ", getResult(G, ff, k, S5, s))
        end
        println(lg, "TopBetweenness : ", getResult(G, ff, k, S6, s))
        println(lg, "TopPageRank : ", getResult(G, ff, k, S7, s))
        println(lg, "TopCloseness : ", getResult(G, ff, k, S8, s))
    end
end

runMain(ARGS[1], ARGS[2])
