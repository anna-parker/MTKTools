
"""
convert scaled reassortment rate ρ to true reassortment rate 
for different simulation type
"""
function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / (2*N)
    elseif simtype == :yule
        return ρ / N
    elseif simtype == :flu
    	return ρ * (n/2)^0.2 / N
    else
        @error "Unrecognized `simtype`."
    end
end

"""
get_trees(no_trees, no_lineages; remove=false, debug=false, c=0.75, ρ = 0.05
returns trees as well as true MCCs
simulate a total of `no_trees` trees using the ARGTools package, 
specifying the `lineage_no` determines the number of nodes in each tree
remove - if internal branches should be removed (i.e. if trees should not be fully resolved, see parameter c)
c - Parameter to describe how resolved trees are
ρ - Reassortment rate scaled to coalescence rate
s - Fraction of leaves that are not sampled at time 0
"""
function get_trees(no_trees, no_lineages; remove=false, c=0.75, ρ = 0.05, N = 10_000, simtype = :flu, s=0.0)
    # Parameters of the ARG simulation
    N = N # pop size
    nmax = no_lineages  # total number of lineages
    if s!=0.0 #if additional parameters should be added
        n0 = 0.25*nmax
        s = s*(nmax-1)*nmax^0.2 /(2*N) # add leaves before time 0 at a rate relative to coalescence
        r = get_r(ρ, 0.75*nmax, N, simtype) # Absolute reassortment rate
    else
        n0 = nmax
        r = get_r(ρ, nmax, N, simtype) # Absolute reassortment rate
    end

    # Simulating the ARG
    arg = ARGTools.SimulateARG.simulate(N, r, n0; s, nmax, K=no_trees, simtype);
    # The trees for the 2 segments
    trees = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);
    if remove
        trees = remove_branches(trees; c=c, N = N)
    end
    return trees, arg
end

function remove_branches(input_trees; c=0.75, N = 10_000)
    trees = [copy(t) for t in input_trees]
    no_trees = length(trees)
    for i in range(1, no_trees)
        tree = trees[i]
        delete_list = String[]
        for node in internals(tree)
            if !node.isroot
                Pr = exp(-node.tau/(c*N))
                if rand() <= Pr
                    push!(delete_list, node.label)
                end
            end
        end
        for node_label in delete_list
            TreeTools.delete!(trees[i], node_label; delete_time=false)
        end
    end
    trees = [convert(TreeTools.Tree{TreeTools.MiscData}, t) for t in trees]
    return trees
end

function get_c(res_rate, rec; n=100, simtype=:flu)
    if simtype == :flu
        if n==50
            c_values = CSV.read("../influenza_c_values_50.txt", DataFrame)
        elseif n==100
            c_values = CSV.read("../influenza_c_values_100.txt", DataFrame)
        end
        try
            keys = c_values.rec_rate 
            if res_rate == 0.3
                values = c_values."res_0.3" 
            elseif res_rate == 0.35
                values = c_values."res_0.35" 
            elseif res_rate == 0.4
                values = c_values."res_0.4" 
            end
            dict_ = Dict(zip(keys, values))
            return dict_[rec]
        catch
            throw(Error("not calculated"))
        end
    elseif simtype == :kingman
        if n==50
            c_values = [0.184, 0.152, 0.1269]
        elseif n==100
            c_values = [0.121, 0.099, 0.079]
        end
        try
            dict_ = Dict(zip([0.3, 0.35, 0.4], c_values))
            return dict_[res_rate] 
        catch
            throw(Error("not calculated"))
        end
    end
end

function get_real_MCCs(no_trees, arg)
    rMCCs = Vector{Vector{String}}[]
    for k in 2:no_trees
        k_iters = Combinatorics.combinations(1:no_trees, k)
        for combination in k_iters
            push!(rMCCs, ARGTools.MCCs_from_arg(arg, combination...));
        end
    end
    return rMCCs
end
