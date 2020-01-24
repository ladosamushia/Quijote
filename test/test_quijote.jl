include("../triplet_count.jl")

using Test

"""
rmax = 3.0
dr = 1.0
xyzw1 = [[0, 0, 0, 1], [1, 0, 0, 1]]
xyzw2 = [[0, 1, 0, 0.5], [0, 0, 1, 1]]
xyzw3 = [[0, 0, 2, 0.5], [10, 0, 0, 1]]
counts = zeros(3, 3, 3)
tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)
counts_test = zeros(3, 3, 3)
counts_test[2,3,1] = 0.5
counts_test[2,3,3] = 0.25
counts_test[1,2,1] = 0.5
counts_test[1,2,3] = 0.25
@test counts == counts_test

@test neighbouring_indeces(1, 2, 3, 3) == [[3, 1, 2], [1, 1, 2], [2, 1, 2], [3, 2, 2], [1, 2, 2], [2, 2, 2], [3, 3, 2], [1, 3, 2], [2, 3, 2], [3, 1, 3], [1, 1, 3], [2, 1, 3], [3, 2, 3], [1, 2, 3], [2, 2, 3], [3, 3, 3], [1, 3, 3], [2, 3, 3], [3, 1, 1], [1, 1, 1], [2, 1, 1], [3, 2, 1], [1, 2, 1], [2, 2, 1], [3, 3, 1], [1, 3, 1], [2, 3, 1]]
"""



