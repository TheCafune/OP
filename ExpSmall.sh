#!/bin/sh
julia -O3 Main.jl Karate 3,6,30
julia -O3 Main.jl Karate 5,6,30
julia -O3 Main.jl Dolphins 3,6,30
julia -O3 Main.jl Dolphins 5,6,30
julia -O3 Main.jl Diseasome 3,6,30
julia -O3 Main.jl Diseasome 5,6,30
julia -O3 Main.jl Netscience 3,6,30
julia -O3 Main.jl Netscience 5,6,30
