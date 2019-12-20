using Test
using DelimitedFiles
include("../triplet_count.jl")

ddd_count_Patchy("test/test_cat_grid_2.txt", "test/test_ddd.txt", 20, 20, 20, 0, 2)
ddr_count_Patchy("test/test_cat_grid_2.txt", "test/test_cat_grid_2.txt", "test/test_ddr.txt", 20, 20, 20, 0, 2)


aa = readdlm("test/test_ddd.txt", comments=true)
bb = readdlm("test/test_ddr.txt", comments=true)
cc = readdlm("test/test_out_grid_2.txt", comments=true)

@test aa == bb
@test aa == cc
