"""
neighbouring_indeces(i, j, k, Nsub)

return indeces of all neighbouring blocks (27 in 3D, including the block itself) for the block (i, j, k).
Nsub is the total number of blocks in each direction.

Function must wrap-around (periodic boundary conditions).
"""
function neighbouring_indeces(i, j, k, Nsub)
    # All neighbours + itself
    i_n = reshape(collect.(Iterators.product(i-1:i+1, j-1:j+1, k-1:k+1)), 27)
    # Periodic boundary conditions
    for index in i_n
        index[index .== 0] .= Nsub
        index[index .== Nsub + 1] .= 1
    end
    return i_n
end

"""
tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)

bin distances between particles in three arrays and incriment histogram in counts.
dr - bin width, rmax - maximum separation.
xyzw - an array of x, y, z, and weight.
"""
function tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)
    for i1 in 1:length(xyzw1)
            for i2 in 1:length(xyzw2)
                r12 = sqrt(sum((xyzw1[i1][1:3] - xyzw2[i2][1:3]).^2))
            if r12 > rmax
                continue
            end
            for i3 in 1:length(xyzw3)
                r13 = sqrt(sum((xyzw1[i1][1:3] - xyzw3[i3][1:3]).^2))
                if r13 > rmax
                    continue
                end
                r23 = sqrt(sum((xyzw2[i2][1:3] - xyzw3[i3][1:3]).^2))
                if r13 > rmax
                    continue
                end
                index = [ceil(Int, r12/dr), ceil(Int, r13/dr), ceil(Int, r23/dr)]
                counts[index[1],index[2],index[3]] += xyzw1[i1][4]*xyzw2[i2][4]*xyzw3[i3][4]
            end
        end
    end
    return nothing
end

"""
cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)

Triple loop over all subcubes and their immediate neighbours
"""
function cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)
    for i1 = 1:Nsub, j1 = 1:Nsub, k1 = 1:Nsub
        i_neighbour = neighbouring_indeces(i1, j1, k1, Nsub)
        xyzw1 = xyzw_cube[i1, j1, k1]
        for ijk2 in length(i_neighbour)
            i2, j2, k2 = i_neighbour[ijk2]
            xyzw2 = xyzw_cube[i2, j2, k2]
            for ijk3 in ijk2:length(i_neighbour)
                i3, j3, k3 = i_neighbour[ijk3]
                xyzw3 = xyzw_cube[i3, j3, k3]
                tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)
            end
        end
    end 
end

function test_triplet_counts()
    # Create random arrays
    Ngal = 10
    Lsurvey = 1000
    xyzw = rand(4, Ngal)
    xyzw[1:3,:] *= Lsurvey

    # Sub-Volumes
    Lsub = 500
    Nsub = ceil(Int, Lsurvey/Lsub)
    println("Nsub ", Nsub)
    # Binning
    rmin = 0
    rmax = 10
    Nbin = 10
    dr = (rmax - rmin)/Nbin

    xyzw_cube = Array{Array{Array{Float64,1}}}(undef, Nsub, Nsub, Nsub)
    i_xyz = ceil.(Int, xyzw[1:3,:]/Lsub)
    for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
        xyzw_cube[i, j, k] = [zeros(4),]
    end
    for i = 1:Ngal
        push!(xyzw_cube[i_xyz[1,i],i_xyz[2,i],i_xyz[3,i]], xyzw[:,i])
    end

    # cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)
end