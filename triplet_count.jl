function neighbouring_indeces(i, j, k, Nsub)
    i_n = reshape(collect.(Iterators.product(i-1:i+1, j-1:j+1, k-1:k+1)), (1, 27))
    for index in i_n
        index[index .== 0] .= Nsub
        index[index .== Nsub + 1] .= 1
    end
    return i_n
end

function test_triplet_counts()
# Create random arrays
Ngal = 1000000
Lsurvey = 1000
xyzw = rand(4, Ngal)
xyzw[1:3,:] *= Lsurvey

# Sub-Volumes
Lsub = 20
Nsub = ceil(Int, Lsurvey/Lsub)

xyzw_cube = Array{Array{Array{Float32}}}(undef, Nsub, Nsub, Nsub)

i_xyz = ceil.(Int, xyzw[1:3,:]/Lsub)

for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
    xyzw_cube[i, j, k] = [[0, 0, 0, 0],]
end

for i = 1:Ngal
    push!(xyzw_cube[i_xyz[i]], xyzw[:,i])
end

#=
# Triplet count
for i1 = 1:Nsub, j1 = 1:Nsub, k1 = 1:Nsub
    xyzw1 = xyzw_cube[i1, j1, k1]
    i_neighbour = neighbouring_indeces(i1, j1, k1, Nsub)
    for ijk2 in Iterators.product(i_neighbour, i_neighbour)
        xyzw2 = xyzw_cube[ijk2[1][1], ijk2[1][2], ijk2[1][3]]
        xyzw3 = xyzw_cube[ijk2[2][1], ijk2[2][2], ijk2[2][3]]
        d12 = sum((xyzw1[1:3] - xyzw2[1:3]).^2)
    end 
end
=#
end

test_triplet_counts()

