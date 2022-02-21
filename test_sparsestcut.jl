include("sparsestcut.jl")

##
mat = matread("Netscience_xy.mat")
A = mat["A"]


n = size(A, 1)
T = Int(ceil(10 * log(n)^2))

min_S, min_expansion = cut_matching(A, T)
@show min_S
@show min_expansion
@show compute_expansion(A, min_S)
##
using MatrixMarket

A = MatrixMarket.mmread("ca-AstroPh.mtx")
n = size(A, 1)
T = Int(ceil(10 * log(n)^2))

min_S, min_expansion = cut_matching(A, T)
@show min_S
@show min_expansion
@show compute_expansion(A, min_S)
##