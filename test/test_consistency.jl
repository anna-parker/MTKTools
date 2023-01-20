using Test
using TreeTools
using TreeKnit
using MTKTools

println("##### testing consistency constraints #####")

function join_sets_slow(input_sets::Vector{Vector{Vector{String}}})
    start_set = input_sets[1]
    for i in 2:length(input_sets)
        joint_sets = Vector{String}[]
        to_be_joint_set = input_sets[i]
        for s1 in start_set
            for s2 in to_be_joint_set
                if !isdisjoint(s1, s2)
                    push!(joint_sets, intersect(s1, s2))
                end
            end
        end
        start_set = joint_sets
    end
    return TreeKnit.sort(start_set; lt=TreeKnit.clt)
end


tree3 = convert(Tree{TreeTools.MiscData}, parse_newick_string("((A,B)i1,((C,D)i2,((E,(F1,F2)i6)i4,G)i5)i3)i7;"))
tree1 = convert(Tree{TreeTools.MiscData}, parse_newick_string("((A,B)j1,(((C,D)j2,(E,(F1,F2)j3)j4)j5,G)j6)j7;"))
tree2 = convert(Tree{TreeTools.MiscData}, parse_newick_string("((G,(A,B)k1)k2,((E,(C,D)k3)k4,(F1,F2)k5)k6)k7;"))
MCC12 = TreeKnit.sort([["A", "B", "E", "F1", "F2"], ["G"], ["C", "D"]], lt=TreeKnit.clt)
MCC13 = TreeKnit.sort([["A", "B", "C", "D", "G"], ["E", "F1", "F2"]], lt=TreeKnit.clt)

constraint = MTKTools.MCC_join_constraint([MCC12, MCC13])

@testset "MCC join" begin
    @test join_sets_slow([MCC12, MCC13]) ==  MTKTools.MCC_join_constraint([MCC12, MCC13])
    @test constraint == [["G"], ["A", "B"], ["C", "D"], ["E", "F1", "F2"]]
end

@testset "mcc_map functions" begin
    mcc_map = TreeKnit.map_mccs(constraint)
    @test mcc_map == Dict("B" => 2, "A" => 2, "C" => 3, "D" => 3, "G" => 1, "E" => 4, "F1" => 4, "F2" => 4)
    TreeKnit.map_mccs!(tree3, constraint)
    inner_node_dict = Dict("i1" => 2, "i2" => 3, "i3" => nothing, "i4" => 4, "i5" => nothing, "i6" => 4, "i7" => nothing)
    for node in nodes(tree3)
        if isleaf(node)
            @test node.data["mcc"] == mcc_map[node.label]
        else
            @test node.data["mcc"] == inner_node_dict[node.label]
        end
    end
end

@testset "map_shared_branches! functions" begin
    MTKTools.map_shared_branches!(constraint, tree3)
    true_labels = Set(["A", "B", "C", "D", "E", "F1", "F2", "i6"])
    for node in nodes(tree3)
        if node.label in true_labels
            @test node.data["shared_branch"] == true
        else
            @test node.data["shared_branch"] == false
        end
    end
end

@testset "map_shared_branches functions" begin
    shared_branches_map_ = MTKTools.map_shared_branches(constraint, tree3)
    true_labels = Set(["A", "B", "C", "D", "E", "F1", "F2", "i6"])
    for node in nodes(tree3)
        if node.label in true_labels
            @test shared_branches_map_[node.label] == true
        else
            @test shared_branches_map_[node.label] == false
        end
    end
end