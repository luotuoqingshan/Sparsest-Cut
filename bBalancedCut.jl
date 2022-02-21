include("sparsestcut.jl")
using Random

function flow_decomposition_with_fake_edges(F::SparseMatrixCSC,
     s::Int64, t::Int64, S::Vector{Int64}, flowvalue::Int64)
    N = size(F, 1) - 2
    S_bar = set_complement(N, S)
    matching_size = div(N, 2) 

    M = zeros(Int, matching_size, 2)

    in_S = zeros(N + 2)
    in_S_bar = zeros(N + 2)
    in_S[S .+ 1] .= 1
    in_S_bar[S_bar .+ 1] .= 1

    Neighs, d = ConstructAdj(F, N + 2) 
    Neighs_ptrs = ones(Int, N + 2)

    visited = zeros(N + 2)

    q = Queue{Pair{Int, Int}}()
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

        v = uNeighs[Neighs_ptrs[u]]
        while F[u, v] <= 0
            Neighs_ptrs[u] += 1
            if Neighs_ptrs[u] > d[u]
                break
            end
            v = uNeighs[Neighs_ptrs[u]]
        end
        if Neighs_ptrs[u] > d[u]
            continue
        end
        if from == s 
            @assert in_S[v] == 1
            from = v
        end
        F[u, v] -= 1
        if v == t 
            edges += 1    
            M[edges, 1] = from 
            M[edges, 2] = u 
            visited[from] = 1
            visited[u] = 1
        else
            enqueue!(q, v => from)
        end
    end
    @assert edges == flowvalue println(string(edges) * " " * string(matching_size))

    unvisited_S = []
    unvisited_S_bar = []
    for u = 2:N + 1
        if visited[u] == 0
            if in_S[u] == 1
                push!(unvisited_S, u)
            else
                push!(unvisited_S_bar, u)
            end
        end
    end
    shuffle(unvisited_S)
    shuffle(unvisited_S_bar)
    for i = 1:matching_size - flowvalue 
        edges += 1
        M[edges, 1] = unvisited_S[i]
        M[edges, 2] = unvisited_S_bar[i]
    end
    return M
end

function flow_match_with_fake_edges(A::SparseMatrixCSC, S::Vector{Int}, beta::Float64)
    N = size(A, 1)
    target_flow = div(N, 2) - ceil(beta * N / log2(N) / log2(N)) 

    lb = 0
    ub = Int(ceil(log2(N)))
    while lb < ub
        mid = div(lb + ub, 2)
        B, s, t = build_graph(A, S, 2^mid)
        F = maxflow(B, s, t)
        if F.flowvalue >= target_flow 
            ub = mid
        else
            lb = mid + 1
        end
    end
    B, s, t = build_graph(A, S, 2^ub)
    F = maxflow(B, s, t)
    M = flow_decomposition_with_fake_edges(F.F, s, t, S, Int(F.flowvalue)) .- 1

    power = max(0, ub - 1)
    B, s, t = build_graph(A, S, 2^power)
    F = maxflow(B, s, t)
    cut_S = source_nodes(F) .- 1
    cut_S = setdiff(cut_S, 0)
    return cut_S, M
end

function pseudo_balanced_cut_matching(A::SparseMatrixCSC, T::Int, b::Float64)
    Expender = expander(Vector{matching}())
    N = size(A, 1) 
    min_expansion = N^2 # store the minimum expansion 
    min_S = [] # store the set inducing the sparsest cut
    for i = 1:N
        S = [i]
        expansion = compute_expansion(A, S) 
        if expansion < min_expansion
            min_expansion = expansion
            min_S = S
        end
    end 
    for _ = 1:T
        r = gen_rand_unitv(N)
        bisec_S = find_bisection(Expender, r)    
        beta = b / 2
        cut_S, M = flow_match_with_fake_edges(A, bisec_S, beta)
        if length(cut_S) > 0 && length(cut_S) < N
            expansion = compute_expansion(A, cut_S)
            if expansion < min_expansion
                min_expansion = expansion
                min_S = cut_S
            end
        end
        push!(Expender.matchings, matching(M))
    end
    return min_S, min_expansion  
end

