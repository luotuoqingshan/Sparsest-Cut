include("bBalancedCut.jl")

##
mat = matread("Netscience_xy.mat")
A = mat["A"]


n = size(A, 1)
T = Int(ceil(10 * log(n)^2))
b = 1/3

min_S, min_expansion = pseudo_balanced_cut_matching(A, T, b)
@show min_S
@show length(min_S) / n
@show min_expansion
@show compute_expansion(A, min_S)

##
using MatrixMarket

A = MatrixMarket.mmread("ca-AstroPh.mtx")
n = size(A, 1)
T = Int(ceil(10 * log(n)^2))

min_S, min_expansion = cut_matching(A, T)
@show min_S
@show length(min_S) / n
@show min_expansion
@show compute_expansion(A, min_S)
##