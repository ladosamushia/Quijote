using Base.Threads
using Base.Iterators
using DelimitedFiles
using StaticArrays
using Printf

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
            elseif xyz3 == xyz1
                i3min = i1
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
                sort!(index)
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
function triplet_counts(Ngal, xyz, w, Lsub, rmin, rmax, Nbin)

    min_xyz = [minimum(xyz[1,:]), minimum(xyz[2,:]), minimum(xyz[3,:])]
    size_xyz = xyz .- min_xyz
    Lsurvey = maximum(size_xyz)
    println("Survey size ", Lsurvey)
    # Number of subcubes is Nsub^3
    Nsub = ceil(Int, Lsurvey/Lsub)
    println("Nsub ", Nsub)
    dr = (rmax - rmin)/Nbin
    # Histogram of triplet counts
    counts = zeros(Float32, nthreads(), Nbin, Nbin, Nbin)
    # Fill in the sub-cubes
    xyz_cube = Array{Array{SVector{3,Float64}}}(undef, Nsub, Nsub, Nsub)
    w_cube = Array{Array{Float64,1}}(undef, Nsub, Nsub, Nsub)
    # Subcube indeces for all galaxies
    i_xyz[i_xyz .== 0] .= 1
    # Fill in the subcubes
    println(w[5], xyz[:,5])
    for i in 1:Ngal
        i1, j1, k1 = i_xyz[:,i]
        if isassigned(xyz_cube, i1, j1, k1)
            push!(xyz_cube[i1, j1, k1], xyz[:,i])
            push!(w_cube[i1, j1, k1], w[i])
        else
            xyz_cube[i1, j1, k1] = [xyz[:,i],]
            w_cube[i1, j1, k1] = [w[i],]
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

function write_counts(ofilename, counts, rmin, rmax, Nbin)
    dr = (rmax - rmin)/Nbin
    ofile = open(ofilename, "w")
    for i in 1:Nbin, j in i:Nbin, k in j:Nbin
        if i + j >= k
            @printf(ofile, "%lf %lf %lf %lf\n", dr*(i - 0.5), dr*(j - 0.5), dr*(k - 0.5), counts[1,i,j,k])
        end
    end
    close(ofile)
    return nothing
end
        


function triplet_count_Patchy(ifilename, ofilename, Lsub, zmin, zmax)
    println(ifilename, " ", Lsub, " ", zmin, " ", zmax)
    patchy = readdlm(ifilename)
    xyz = patchy[:,1:3]
    w = patchy[:,4]
    for i in 1:length(w)
        if w[i] != 0
            w[i] = 1
        end
    end
    xyz = transpose(xyz)
    w = transpose(w)
    Ngal = size(w)[2]
    println("Ngal ", Ngal)
    counts = triplet_counts(Ngal, xyz, w, Lsub, 0, 20, 20)
    write_counts(ofilename, counts, 0, 20, 20)
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end
