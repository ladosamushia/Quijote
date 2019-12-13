using Base.Threads
using Base.Iterators
using DelimitedFiles

"""
unique_triplet_indeces()

This function computes unique triplets I could go to from any point.
"""
function unique_triplet_indeces()

    i_n = Array{Int,1}[]
    # Only move forward
    for di in -1:1, dj in -1:1
        push!(i_n, [di, dj, 1])
        if dj >= di && dj > -1
            push!(i_n, [di, dj, 0])
        end
    end
    # Now triplets
    j_n = []
    for i in 1:length(i_n), j in i:length(i_n)
        # Are the two neighbours themselves neighbours?
        if maximum(abs.(i_n[i] - i_n[j])) <= 1
            push!(j_n, [i_n[i], i_n[j]])
        end
    end
    return j_n

end

"""
tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)

bin distances between particles in three arrays and incriment histogram in counts.

dr - bin width, rmax - maximum separation.
xyzw - an array of x, y, z, and weight.
"""

function tri_bin(xyz1, xyz2, xyz3, w1, w2, w3, dr, rmax, counts)
    
    for i1 in 1:length(xyz1)
        if xyz1 == xyz2
            i2min = i1
        else
            i2min = 1
        end
        for i2 in i2min:length(xyz2)
            r12 = sqrt(sum((xyz1[i1] - xyz2[i2]).^2))
            if r12 >= rmax || r12 <= 1e-3
                continue
            end
            if xyz2 == xyz3
                i3min = i2
            else
                i3min = 1
            end
            for i3 in i3min:length(xyz3)
                r13 = sqrt(sum((xyz1[i1] - xyz3[i3]).^2))
                if r13 >= rmax || r13 <= 1e-3
                    continue
                end
                r23 = sqrt(sum((xyz2[i2] - xyz3[i3]).^2))
                if r23 >= rmax || r23 <= 1e-3
                    continue
                end
                index = [ceil(Int, r12/dr), ceil(Int, r13/dr), ceil(Int, r23/dr)]
                counts[threadid(),index[1],index[2],index[3]] += w1[i1]*w2[i2]*w3[i3]
                
            end
        end
    end
    return nothing

end

"""
cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)

Triple loop over all subcubes and their immediate neighbours
"""

function cube_triplets(xyz_cube, w_cube, Nsub, dr, rmax, counts)

    println("start cube_triplets")
    j_n = unique_triplet_indeces()
    # All subcubes
    @threads for i1 in 1:Nsub
        for j1 in 1:Nsub
            for k1 in 1:Nsub
                # Skip empty subcubes
                if isassigned(xyz_cube, i1, j1, k1)
                    xyz1 = xyz_cube[i1, j1, k1]
                    w1 = w_cube[i1, j1, k1]
                else
                    continue
                end
                # All unique pairs of neighbours
                for index23 in j_n
                    i2, j2, k2 = index23[1] + [i1, j1, k1]
                    i3, j3, k3 = index23[2] + [i1, j1, k1]
                    all_indeces = [i2, j2, k2, i3, j3, k3]
                    # Make sure we are inside the cube
                    if maximum(all_indeces) <=Nsub && minimum(all_indeces) >= 1
                        # Skip empty subcubes
                        if isassigned(xyz_cube, i2, j2, k2) && isassigned(xyz_cube, i3, j3, k3)
                            xyz2 = xyz_cube[i2, j2, k2]
                            w2 = w_cube[i2, j2, k2]
                            xyz3 = xyz_cube[i3, j3, k3]
                            w3 = w_cube[i3, j3, k3]
                            tri_bin(xyz1, xyz2, xyz3, w1, w2, w3, dr, rmax, counts)
                        else
                            continue
                        end
                    end
                end
            end
        end
    end

end

"""
triplet_counts(Ngal, Lsurvey, xyzw, Lsub, rmin, rmax, Nbin)

Count all triplets.


"""

function triplet_counts(Ngal, Lsurvey, xyzw, Lsub, rmin, rmax, Nbin)

    # Number of subcubes is Nsub^3
    Nsub = ceil(Int, Lsurvey/Lsub)
    dr = (rmax - rmin)/Nbin
    # Histogram of triplet counts
    counts = zeros(Float32, nthreads(), Nbin, Nbin, Nbin)
    # Fill in the sub-cubes
    xyz_cube = Array{Array{Array{Float64,1}}}(undef, Nsub, Nsub, Nsub)
    w_cube = Array{Array{Float64,1}}(undef, Nsub, Nsub, Nsub)
    min_xyz = [minimum(xyzw[1,:]), minimum(xyzw[2,:]), minimum(xyzw[3,:])]
    # Subcube indeces for all galaxies
    i_xyz = ceil.(Int, (xyzw[1:3,:] .- min_xyz)/Lsub)
    i_xyz[i_xyz .== 0] .= 1
    # Fill in the subcubes
    for i in 1:Ngal
        i1, j1, k1 = i_xyz[:,i]
        if isassigned(xyz_cube, i1, j1, k1)
            push!(xyz_cube[i1, j1, k1], xyzw[1:3,i])
            push!(w_cube[i1, j1, k1], xyzw[4,i])
        else
            xyz_cube[i1, j1, k1] = [xyzw[1:3,i],]
            w_cube[i1, j1, k1] = [xyzw[4,i],]
        end
    end
    # Count triplets
    println("cube_triplets")
    cube_triplets(xyz_cube, w_cube, Nsub, dr, rmax, counts)
    # Reduce counts across threads
    counts = sum(counts, dims=1)
    println("sum", sum(counts))
    return counts

end

function test_triplet_counts()
    xyzw = rand(4, 1000000)
    xyzw[1:3,:] *= 1000
    triplet_counts(1000000, 1000, xyzw, 20, 0, 20, 20)
    return nothing
end

function test_on_Patchy(stride, Lsub)
    patchy = readdlm("Patchy-Mocks-DR12NGC-COMPSAM_V6C_0001_xyzw.dat")
    xyzw = patchy[1:stride:end,1:4]
    xyzw[:,4] .= 1
    xyzw = transpose(xyzw)
    Ngal = size(xyzw)[2]
    println("Ngal ", Ngal)
    Lsurvey = 3300.0
    triplet_counts(Ngal, Lsurvey, xyzw, Lsub, 0, 20, 20)
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end

test_on_Patchy(1,20)
