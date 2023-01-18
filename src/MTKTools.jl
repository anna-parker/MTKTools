module MTKTools

using TreeTools
using TreeKnit
using Combinatorics
using DataFrames
using CSV
using ARGTools

include("simulations.jl")
include("consistency.jl")
export map_shared_branches, map_shared_branches!, is_MCC_subset, MCC_join_constraint
include("measures.jl")
export get_trees, get_c, get_real_MCCs

end