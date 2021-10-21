#!/bin/sh
julia -O3 Main.jl USgrid 5,100,5000
julia -O3 Main.jl USgrid 10,100,5000
julia -O3 Main.jl GrQc 5,100,5000
julia -O3 Main.jl GrQc 10,100,5000
julia -O3 Main.jl HepPh 5,100,5000
julia -O3 Main.jl HepPh 10,100,5000
julia -O3 Main.jl Gplus 5,100,5000
julia -O3 Main.jl Gplus 10,100,5000
