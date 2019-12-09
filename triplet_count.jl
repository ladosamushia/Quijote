using Base.Threads

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
    for i in 1:27
        for j in 1:3
            if i_n[i][j] == 0
                i_n[i][j] = Nsub
            elseif i_n[i][j] == Nsub + 1
                i_n[i][j] = 1
            end
        end
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
            if r12 >= rmax || r12 == 0
                continue
            end
            for i3 in 1:length(xyzw3)
                r13 = sqrt(sum((xyzw1[i1][1:3] - xyzw3[i3][1:3]).^2))
                if r13 >= rmax || r13 == 0
                    continue
                end
                r23 = sqrt(sum((xyzw2[i2][1:3] - xyzw3[i3][1:3]).^2))
                if r23 >= rmax || r23 == 0
                    continue
                end
                index = [ceil(Int, r12/dr), ceil(Int, r13/dr), ceil(Int, r23/dr)]
                counts[threadid(),index[1],index[2],index[3]] += xyzw1[i1][4]*xyzw2[i2][4]*xyzw3[i3][4]
                
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

   @threads for i1 = 1:Nsub, j1 = 1:Nsub, k1 = 1:Nsub
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

function triplet_counts(Ngal, Lsurvey, xyzw, Lsub, rmin, rmax, Nbin)

    Nsub = ceil(Int, Lsurvey/Lsub)
    println(Nsub)
    dr = (rmax - rmin)/Nbin
    println(dr)
    println(nthreads())
    counts = zeros(Float32, nthreads(), Nbin, Nbin, Nbin)
    println(typeof(counts))
    # Fill in the sub-cubes
    xyzw_cube = Array{Array{Array{Float64,1}}}(undef, Nsub, Nsub, Nsub)
    i_xyz = ceil.(Int, xyzw[1:3,:]/Lsub)
    for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
        xyzw_cube[i, j, k] = [zeros(4),]
    end
    for i in 1:Ngal
        push!(xyzw_cube[i_xyz[1,i],i_xyz[2,i],i_xyz[3,i]], xyzw[:,i])
    end
    # Count triplets
    println("cube_triplets")
    cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)
    counts = sum(counts, dims=1)
    println(counts)
    println(typeof(counts))
    return counts

end

xyzw = rand(4, 10)
xyzw[1:3,:] *= 1000
triplet_counts(10, 1000, xyzw, 20, 0, 20, 20)

"""
function test_triplet_counts()
    # Create random arrays
    Ngal = 1000000
    Lsurvey = 1000
    xyzw = rand(4, Ngal)
    xyzw[1:3,:] *= Lsurvey
    # Sub-Volumes
    Lsub = 20
    Nsub = ceil(Int, Lsurvey/Lsub)
    println("Nsub ", Nsub)
    # Binning
    rmin = 0
    rmax = 20
    Nbin = 20
    dr = (rmax - rmin)/Nbin
    counts = zeros(Float32, Nbin, Nbin, Nbin)
    # Fill in the sub-cubes
    xyzw_cube = Array{Array{Array{Float64,1}}}(undef, Nsub, Nsub, Nsub)
    i_xyz = ceil.(Int, xyzw[1:3,:]/Lsub)
    for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
        xyzw_cube[i, j, k] = [zeros(4),]
    end
    for i = 1:Ngal
        push!(xyzw_cube[i_xyz[1,i],i_xyz[2,i],i_xyz[3,i]], xyzw[:,i])
    end
    # Count triplets
    cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)
end
"""
