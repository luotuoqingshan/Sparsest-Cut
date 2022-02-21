include("maxflow.jl")

using MatrixNetworks
using SparseArrays
using LinearAlgebra
using DataStructures
using MAT

# a list of edges which forms a matching
# size of matrix edges is (matching_size, 2)
mutable struct matching
    edges::Matrix{Int64}
end

# a list of matchings whose union formes the expander
# each matching is a matrix of size [div(n, 2), 2]
mutable struct expander
    matchings::Vector{matching} end

"""
Generate an unit vector in R^n which is orthogonal to all ones vector
n: the dimension of the vector

Return r, which is a vector satisfying r' * ones(n) = zeros(n) 
"""
function gen_rand_unitv(n::Int)
    # Generate one vector which is element-wise Gaussian 
    r = randn((n))

    # remove the component which is on the direction of ones(n)
    vones = ones(n)
    r = r - r' * vones * vones / n

    # return the normalized vector
    return r / norm(r)
end

"""
Function for cut player to choose a bisection based on previous matchings
returned by matching player.

H = collections of matchings returned by matching player in previous rounds
r = a random unit vector which is orthogonal to all-ones vector 

Return:
S = a set of vertices which induce the bisection 
"""
function find_bisection(H::expander, r::Vector{Float64})
    n = length(r)
    newr = zeros(n)

    # Here we Compute M_i(M_{i-1}(...(M_2(M_1(r)))))
    # where M(r) means projecting distribution r onto matching M, i.e. 
    # do a one-step lazy random walk on M starting from distribution r 
    for i = 1:length(H.matchings)
        M = H.matchings[i] 
        for j = 1:size(M.edges, 1)
            u = M.edges[j, 1]
            v = M.edges[j, 2]

            # for a one-step lazy random walk, w.p. 1/2 it will stay 
            #                                  w.p. 1/2 it will traverse to the neighbor
            newr[u] = (r[u] + r[v]) / 2
            newr[v] = newr[v]
        end
        r = newr
    end

    # sort the vertices according to r
    # intuitively, if in \cup_i M_i two vertices are well-connected(a lot of paths between),
    # then their scores will be mixed well during the projections.
    # Therefore by splitting them according to score, it's more possible to find 
    # some bisection which induces a sparse cut
    id = sortperm(r)
    S = id[1:div(n, 2)]
    return S
end

"""
Given a set S which is a subset of [N], return [N] - S
"""
function set_complement(N::Int, S::Vector{Int})
    V = collect(1:N)
    S_bar = setdiff(V, S)
    return S_bar
end

"""
Given the adjacency matrix A of a graph G and a seed set S,
return the expansion of S w.r.t. G.

Here Expansion(S, G) := E(S, V(G) - S) / min(|S|, |V(G) - S|)
"""
function compute_expansion(A::SparseMatrixCSC, S::Vector{Int64}) 
    N = size(A, 1)
    S_bar = set_complement(N, S) 

    # find the non-zeros entries of Sparse Adjacency Matrix A indexed by [S, S_bar]
    # which is the number of edges of the cut (S, S_bar)
    _, _, E = findnz(A[S, S_bar])
    return sum(E) / min(length(S), length(S_bar))
end

"""
Given the bisection (S, S_bar) of G, we build a graph G' where we assign capacity 
Cap to each edge in G and add one super source s, one super sink t. Further we 
add one undirected edge from s to each vertex v in S with capacity 1, one undirected 
edge from each vertex v in S_bar to t with capacity 1 as well.

A = N by N adjacnecy matrix of undirected graph G
S = a seed set of vertices, which is a subset of vertices of G
Cap = the capacity we assign to each edge of G when building G'

Return:
C = the capacity matrix of G'
s = super sourse's index
t = super sink's index

"""
function build_graph(A::SparseMatrixCSC, S::Vector{Int64}, Cap::Int64)
    N = size(A, 1)
    S_bar = set_complement(N, S)

    # Here for simplicity we let s = 1, t = N + 2 and shift the index of each vertex
    # of G by 1.
    s = 1
    t = N + 2
    svec = zeros(N)
    svec[S] .= 1
    tvec = zeros(N)
    tvec[S_bar] .= 1
    C = [spzeros(1, 1) sparse(svec') spzeros(1, 1);
        sparse(svec) Cap * A sparse(tvec);
        spzeros(1, 1) sparse(tvec') spzeros(1, 1)]
    return C, s, t
end

"""
Given a max s-t flow and a bisection (S, S_bar) of original graph G, we decompose
the s-t flow into a set of s-t paths. Further by picking the first vertex belonging
to S and the last vertex belonging to S_bar, we map each s-t path to a path between
one vertex from S and one vertex from S_bar. Thus we get a matching between 
S and S_bar which is routable in G.

F = Flow matrix of the max s-t flow. 
    F[u, v] > 0 means there exists F[u, v] unit flow flowing from u to v
    F[u, v] < 0 means there exists -F[u, v] unit flow flowing from v to u

s = source's index
t = sink's index
S = a set of vertices which induces the bisection

Return:
M =  a matching between S and S_bar which is routable in G 

Notes: For convenience, we let s = 1 and t = N + 2 when building auxiliary graphs
for computing matching. So here it's always true that s = 1, t = N + 2.
"""
function flow_decomposition(F::SparseMatrixCSC, s::Int64, t::Int64, S::Vector{Int64})
    N = size(F, 1) - 2
    S_bar = set_complement(N, S)

    # target matching size = N // 2
    matching_size = div(N, 2)

    # store the result matching
    M = zeros(Int, matching_size, 2)

    # indicator variables for S and S_bar
    in_S = zeros(N + 2)
    in_S_bar = zeros(N + 2)
    in_S[S .+ 1] .= 1
    in_S_bar[S_bar .+ 1] .= 1

    # Compute Neighbors and degree for each vertex based on the flow matrix
    # Notice here neighbors include both in and out neighbors
    Neighs, d = ConstructAdj(F, N + 2) 

    # pointers for Neighbor lists used for BFS
    Neighs_ptrs = ones(Int, N + 2)

    # queue used for BFS, each item is a pair, which is 
    # (starting vertex, ending vertex) of a path 
    q = Queue{Pair{Int, Int}}()

    # initialize each path with (s, s)
    for _ = 1:length(S)
        enqueue!(q, s => s)
    end

    edges = 0
    while !isempty(q)
        x = dequeue!(q)
        u = x.first
        from = x.second
        uNeighs = Neighs[u]
        @assert Neighs_ptrs[u] <= d[u] 

        # find the first outgoing neighbor with positive flow remaining
        v = uNeighs[Neighs_ptrs[u]]
        while F[u, v] <= 0
            Neighs_ptrs[u] += 1
            @assert Neighs_ptrs[u] <= d[u]
            v = uNeighs[Neighs_ptrs[u]]
        end

        # since finally we want to match vertices in S with vertices in S_bar
        # we change the starting vertex at the first time we meet any vertex in S
        if from == s 
            @assert in_S[v] == 1
            from = v
        end

        # use one unit flow
        F[u, v] -= 1

        # if we reach the sink, it means we have already find one s-t path,
        # store it in M
        if v == t 
            edges += 1    
            M[edges, 1] = from 
            M[edges, 2] = u 
        else
        # if not, then we push the temporary path into the queue again
            enqueue!(q, v => from)
        end
    end
    @assert edges == matching_size println(string(edges) * " " * string(matching_size))
    return M
end

"""
Given the bisection (S, S_bar), binary search capacity c such that 
the matching between (S, S_bar) can be routed in G when we assign each
edge with capacity c and can't be routed when the capacity is c/2.

The max s-t flow we get when the capacity of each edge is assigned to 
c can give us one G-routable matching between (S, S_bar).

The min s-t cut according to the max s-t flow we get when the capacity
of each edge is assigned to c/2 can give us one sparse cut with expansion
at most 2/c.

A = adjacency matrix of G
S = a set of vertices which induces the bisection

Return:
cut_S = a set of vertices which induces sparse cut with expansion at most 2/c 
M = a matching between (S, S_bar) which is routable in G when capacity is c
"""
function flow_match(A::SparseMatrixCSC, S::Vector{Int})
    N = size(A, 1)
    lb = 0
    ub = Int(ceil(log2(N)))

    # binary search the smallest capacity (within factor 2) c such that
    # a matching between (S, S_bar) can be routed in G 
    while lb < ub
        mid = div(lb + ub, 2)
        B, s, t = build_graph(A, S, 2^mid)
        F = maxflow(B, s, t)
        if F.flowvalue >= div(N, 2) 
            ub = mid
        else
            lb = mid + 1
        end
    end
    # compute the matching
    # Here we need to compute the flow when the capacity of each edge is c
    B, s, t = build_graph(A, S, 2^ub)
    F = maxflow(B, s, t)
    M = flow_decomposition(F.F, s, t, S) .- 1

    # compute the sparse cut
    # Here we need to compute the flow when the capacity of each edge is c/2
    power = max(0, ub - 1)
    B, s, t = build_graph(A, S, 2^power)
    F = maxflow(B, s, t)

    # Let S be the set of vertices reachable from sourse in the residual graph
    # which induces the desired min cut
    cut_S = source_nodes(F) .- 1
    cut_S = setdiff(cut_S, 0)
    return cut_S, M
end

"""
Cut-matching game:
In round i,
Cut Player: Based on the matching M_1, M_2, ..., M_{i-1} returned by matching player
in the first (i-1) rounds. Find a new bisection (S, S_bar) for this round.

Matching Player: Given the bisection (S, S_bar) chosen by the cut player, find
a matching between S and S_bar, which is routable in G with congestion c_i. 
Meanwhile output a sparse cut whose expansion is at most c_i/2. 

After T = Omega(log2(N)^2) rounds, we find one expander which is routable 
in G with congestion (T * max_i c_i), and one sparse cut with expansion at most
(2 / max_i c_i). This certifies that phi(G) (expansion of G) is at most 1/T *   
phi(S).

We can show that the union of M_1, M_2, ..., M_T is an expander with high probability
for some T = O(log2(N)^2)

A = adjacency matrix of G

"""
function cut_matching(A::SparseMatrixCSC, T::Int)
    # store the set of matchings which forms an expander
    Expender = expander(Vector{matching}())
    N = size(A, 1) 

    # store the minimum expansion 
    min_expansion = N^2 

    # store the set inducing the sparsest cut
    min_S = [] 

    # initialize the min_S, min_expansion with S = {v}, i.e. single vertex
    for i = 1:N
        S = [i]
        expansion = compute_expansion(A, S) 
        if expansion < min_expansion
            min_expansion = expansion
            min_S = S
        end
    end 

    # run cut-matching game for T rounds 
    for _ = 1:T
        # Cut Player: give the bisection
        r = gen_rand_unitv(N)
        bisec_S = find_bisection(Expender, r)    

        # Matching Player: give the matching and sparse cut
        cut_S, M = flow_match(A, bisec_S)

        # update the min_expansion and min_S with new sparse cut
        if length(cut_S) > 0 && length(cut_S) < N
            expansion = compute_expansion(A, cut_S)
            if expansion < min_expansion
                min_expansion = expansion
                min_S = cut_S
            end
        end

        #store the new matching into expander
        push!(Expender.matchings, matching(M))
    end
    return min_S, min_expansion  
end

