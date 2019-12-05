function neighbouring_indeces(i, j, k, Nsub)
    i_n = reshape(collect.(Iterators.product(i-1:i+1, j-1:j+1, k-1:k+1)), (1, 27))
    for index in i_n
        index[index .== 0] .= Nsub
        index[index .== Nsub + 1] .= 1
    end
    return i_n
end

function tri_bin(x1, y1, z1, x2, y2, z2, x3, y3, z3, dr, Nbin, rmax)
    r12 = sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2)
    r13 = sqrt((x1 - x3).^2 + (y1 - y3).^2 + (z1 - z3).^2)
    r23 = sqrt((x2 - x3).^2 + (y2 - y3).^2 + (z2 - z3).^2)
end
function test_triplet_counts()
# Create random arrays
Ngal = 10000
Lsurvey = 1000
xyzw = rand(4, Ngal)
xyzw[1:3,:] *= Lsurvey

# Sub-Volumes
Lsub = 250
Nsub = ceil(Int, Lsurvey/Lsub)

x_cube = Array{Array{Float32}}(undef, Nsub, Nsub, Nsub)
y_cube = Array{Array{Float32}}(undef, Nsub, Nsub, Nsub)
z_cube = Array{Array{Float32}}(undef, Nsub, Nsub, Nsub)
w_cube = Array{Array{Float32}}(undef, Nsub, Nsub, Nsub)

i_xyz = ceil.(Int, xyzw[1:3,:]/Lsub)

for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
    x_cube[i, j, k] = [0,]
    y_cube[i, j, k] = [0,]
    z_cube[i, j, k] = [0,]
    w_cube[i, j, k] = [0,]
end

for i = 1:Ngal
    push!(x_cube[i_xyz[i]], xyzw[1,i])
    push!(y_cube[i_xyz[i]], xyzw[2,i])
    push!(z_cube[i_xyz[i]], xyzw[3,i])
    push!(w_cube[i_xyz[i]], xyzw[4,i])
end

# Triplet count
for i1 = 1:Nsub, j1 = 1:Nsub, k1 = 1:Nsub
    x1 = x_cube[i1, j1, k1]
    y1 = y_cube[i1, j1, k1]
    z1 = z_cube[i1, j1, k1]
    w1 = w_cube[i1, j1, k1]
    i_neighbour = neighbouring_indeces(i1, j1, k1, Nsub)
    for ijk2 in Iterators.product(i_neighbour, i_neighbour)
        x2 = x_cube[ijk2[1][1], ijk2[1][2], ijk2[1][3]]
        y2 = y_cube[ijk2[1][1], ijk2[1][2], ijk2[1][3]]
        z2 = z_cube[ijk2[1][1], ijk2[1][2], ijk2[1][3]]
        w2 = w_cube[ijk2[1][1], ijk2[1][2], ijk2[1][3]]
        x3 = x_cube[ijk2[2][1], ijk2[2][2], ijk2[2][3]]
        y3 = y_cube[ijk2[2][1], ijk2[2][2], ijk2[2][3]]
        z3 = z_cube[ijk2[2][1], ijk2[2][2], ijk2[2][3]]
        w3 = w_cube[ijk2[2][1], ijk2[2][2], ijk2[2][3]]
        for ii1 = 1:length(x1)
            dx12 = x1[ii1] .- x2
            dy12 = y1[ii1] .- y2
            dz12 = z1[ii1] .- z2
            dx13 = x1[ii1] .- x3
            dy13 = y1[ii1] .- y3
            dz13 = z1[ii1] .- z3
            r12 = sqrt.(dx12.^2 + dy12.^2 + dz12.^2)
            r23 = sqrt.(dx13.^2 + dy13.^2 + dz13.^2)
        end
    end 
end

end

test_triplet_counts()

