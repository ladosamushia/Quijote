using Base.Threads
using Base.Iterators
using DelimitedFiles
using StaticArrays
using Printf

function unique_neighbouring_cubes!(i_n)
    # Only move forward
    counter = 1
    for di in -1:1, dj in -1:1
        i_n[:,counter] = [di, dj, 1]
        counter = counter + 1
        if dj >= di && dj > -1
            i_n[:,counter] = [di, dj, 0]
            counter += 1
        end
    end
    return nothing
end

"""
unique_triplet_indeces()

This function computes unique triplets I could go to from any point.
"""
function unique_triplet_indeces()
    i_n = zeros(Int, 3, 14)
    unique_neighbouring_cubes!(i_n)
    j_n = zeros(Int, 2, 3, 71)
    # Now triplets
    counter = 1
    for i in 1:14, j in i:14
        # Are the two neighbours themselves neighbours?
        if maximum(abs.(i_n[:,i] - i_n[:,j])) <= 1
            j_n[1,:,counter] = i_n[:,i]
            j_n[2,:,counter] = i_n[:,j]
            counter += 1
        end
    end
    return j_n
end

function double_bin(xyz1, xyz2, w1, w2, dr, rmax, counts)
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
            index = ceil(Int, r12/dr)
            counts[threadid(),index] += w1[i1]*w2[i2]
        end
    end
    return nothing
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
                i12 = ceil(Int, r12/dr)
                i13 = ceil(Int, r13/dr)
                i23 = ceil(Int, r23/dr)
                # Sorting
                imin = min(i12, i13, i23)
                imax = max(i12, i13, i23)
                if i12 <= max(i13, i23) && i12 >= min(i13, i23)
                    imid = i12
                elseif i13 <= max(i12, i23) && i13 >= min(i12, i23)
                    imid = i13
                else
                    imid = i23
                end
                counts[threadid(),imin,imid,imax] += w1[i1]*w2[i2]*w3[i3]
                
            end
        end
    end
    return nothing
end

function cube_pairs(xyz_cube_1, xyz_cube_2, w_cube_1, w_cube_2, Nsub, dr, rmax, counts)
    println("start cube_pairs")
    i_n = zeros(Int, 3, 14)
    unique_neighbouring_cubes!(i_n)
    # All subcubes
    @threads for i1 in 1:Nsub
        for j1 in 1:Nsub
            for k1 in 1:Nsub
                # Skip empty subcubes
                if isassigned(xyz_cube_1, i1, j1, k1)
                    xyz1 = xyz_cube_1[i1, j1, k1]
                    w1 = w_cube_1[i1, j1, k1]
                else
                    continue
                end
                # All unique pairs of neighbours
                for cc in 1:14
                    i2 = j_n[1,1,cc] + i1
                    j2 = j_n[1,2,cc] + j1
                    k2 = j_n[1,3,cc] + k1
                    # Make sure we are inside the cube
                    if max(i2, j2, k2, i3, j3, k3) <=Nsub && min(i2, j2, k2, i3, j3, k3) >= 1
                        # Skip empty subcubes
                        if isassigned(xyz_cube_2, i2, j2, k2)
                            xyz2 = xyz_cube_2[i2, j2, k2]
                            w2 = w_cube_2[i2, j2, k2]
                            double_bin(xyz1, xyz2, w1, w2, dr, rmax, counts)
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
cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)

Triple loop over all subcubes and their immediate neighbours
"""
function cube_triplets(xyz_cube_12, xyz_cube_3, w_cube_12, w_cube_3, Nsub, dr, rmax, counts)
    println("start cube_triplets")
    j_n = unique_triplet_indeces()
    # All subcubes
    @threads for i1 in 1:Nsub
        for j1 in 1:Nsub
            for k1 in 1:Nsub
                # Skip empty subcubes
                if isassigned(xyz_cube_12, i1, j1, k1)
                    xyz1 = xyz_cube_12[i1, j1, k1]
                    w1 = w_cube_12[i1, j1, k1]
                else
                    continue
                end
                # All unique pairs of neighbours
                for cc in 1:71
                    i2 = j_n[1,1,cc] + i1
                    j2 = j_n[1,2,cc] + j1
                    k2 = j_n[1,3,cc] + k1
                    i3 = j_n[2,1,cc] + i1
                    j3 = j_n[2,2,cc] + j1
                    k3 = j_n[2,3,cc] + k1
                    # Make sure we are inside the cube
                    if max(i2, j2, k2, i3, j3, k3) <=Nsub && min(i2, j2, k2, i3, j3, k3) >= 1
                        # Skip empty subcubes
                        if isassigned(xyz_cube_12, i2, j2, k2) && isassigned(xyz_cube_3, i3, j3, k3)
                            xyz2 = xyz_cube_12[i2, j2, k2]
                            w2 = w_cube_12[i2, j2, k2]
                            xyz3 = xyz_cube_3[i3, j3, k3]
                            w3 = w_cube_3[i3, j3, k3]
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
function triplet_counts(xyz_12, xyz_3, w_12, w_3, Lsub, rmin, rmax, Nbin)
    Ngal12 = size(w_12)[2]
    Ngal3 = size(w_3)[2]
    min_xyz_12 = [minimum(xyz_12[1,:]), minimum(xyz_12[2,:]), minimum(xyz_12[3,:])]
    min_xyz_3 = [minimum(xyz_3[1,:]), minimum(xyz_3[2,:]), minimum(xyz_3[3,:])]
    min_xyz = [min(min_xyz_12[1], min_xyz_3[1]), min(min_xyz_12[2], min_xyz_3[2]), min(min_xyz_12[3], min_xyz_3[3])]
    max_xyz_12 = [maximum(xyz_12[1,:]), maximum(xyz_12[2,:]), maximum(xyz_12[3,:])]
    max_xyz_3 = [maximum(xyz_3[1,:]), maximum(xyz_3[2,:]), maximum(xyz_3[3,:])]
    max_xyz = [max(max_xyz_12[1], max_xyz_3[1]), max(max_xyz_12[2], max_xyz_3[2]), max(max_xyz_12[3], max_xyz_3[3])]
    Lsurvey = maximum(max_xyz - min_xyz)
    println("Survey size ", Lsurvey)
    # Number of subcubes is Nsub^3
    Nsub = ceil(Int, Lsurvey/Lsub)
    println("Nsub ", Nsub)
    dr = (rmax - rmin)/Nbin
    # Histogram of triplet counts
    counts = zeros(Float32, nthreads(), Nbin, Nbin, Nbin)
    # Fill in the sub-cubes
    xyz_cube_12 = Array{Array{SVector{3,Float64}}}(undef, Nsub, Nsub, Nsub)
    w_cube_12 = Array{Array{Float64,1}}(undef, Nsub, Nsub, Nsub)
    if xyz_12 != xyz_3
        xyz_cube_3 = Array{Array{SVector{3,Float64}}}(undef, Nsub, Nsub, Nsub)
        w_cube_3 = Array{Array{Float64,1}}(undef, Nsub, Nsub, Nsub)
    else
        xyz_cube_3 = xyz_cube_12
        w_cube_3 = w_cube_12
    end
    # Subcube indeces for all galaxies
    println(size(xyz_12))
    println(size(min_xyz))
    i_xyz = ceil.(Int, (xyz_12 .- min_xyz)/Lsub)
    i_xyz[i_xyz .== 0] .= 1
    # Fill in the subcubes
    for i in 1:Ngal12
        i1, j1, k1 = i_xyz[:,i]
        if isassigned(xyz_cube_12, i1, j1, k1)
            push!(xyz_cube_12[i1, j1, k1], xyz_12[:,i])
            push!(w_cube_12[i1, j1, k1], w_12[i])
        else
            xyz_cube_12[i1, j1, k1] = [xyz_12[:,i],]
            w_cube_12[i1, j1, k1] = [w_12[i],]
        end
    end
    if xyz_12 != xyz_3
        i_xyz = ceil.(Int, (xyz_3 .- min_xyz)/Lsub)
        i_xyz[i_xyz .== 0] .= 1
        for i in 1:Ngal3
            i1, j1, k1 = i_xyz[:,i]
            if isassigned(xyz_cube_3, i1, j1, k1)
                push!(xyz_cube_3[i1, j1, k1], xyz_3[:,i])
                push!(w_cube_3[i1, j1, k1], w_3[i])
            else
                xyz_cube_3[i1, j1, k1] = [xyz_3[:,i],]
                w_cube_3[i1, j1, k1] = [w_3[i],]
            end
        end
    end

    # Count triplets
    println("cube_triplets")
    cube_triplets(xyz_cube_12, xyz_cube_3, w_cube_12, w_cube_3, Nsub, dr, rmax, counts)
    # Reduce counts across threads
    counts = sum(counts, dims=1)
    println("sum", sum(counts))
    return counts
end

function write_counts(ofilename, counts, rmin, rmax, Nbin, Nwgal)
    dr = (rmax - rmin)/Nbin
    ofile = open(ofilename, "w")
    @printf(ofile, "# Weighted number of galaxies: %.1lf\n", Nwgal)
    for i in 1:Nbin, j in i:Nbin, k in j:Nbin
        if i + j >= k
            @printf(ofile, "%.1lf %.1lf %.1lf %.1lf\n", dr*(i - 0.5), dr*(j - 0.5), dr*(k - 0.5), counts[1,i,j,k])
        end
    end
    close(ofile)
    return nothing
end

function ddd_count_Patchy(ifilename, ofilename, Lsub, rmax, Nbin, zmin, zmax)
    println(ifilename, " ", Lsub, " ", zmin, " ", zmax)
    patchy = readdlm(ifilename)
    red = patchy[:,5]
    in_shell = (red .< zmax) .& (red .> zmin)
    xyz = patchy[in_shell,1:3]
    w = patchy[in_shell,4]
    for i in 1:length(w)
        if w[i] != 0
            w[i] = 1
        end
    end
    Nwgal = sum(w)
    xyz = transpose(xyz)
    w = transpose(w)
    counts = triplet_counts(xyz, xyz, w, w, Lsub, 0, rmax, Nbin)
    write_counts(ofilename, counts, 0, 20, 20, Nwgal)
    return nothing
end

function ddr_count_Patchy(ifilename_d, ifilename_r, ofilename, Lsub, rmax, Nbin, zmin, zmax)
    patchy = readdlm(ifilename_d)
    red = patchy[:,5]
    in_shell = (red .< zmax) .& (red .> zmin)
    xyz_12 = patchy[in_shell,1:3]
    w_12 = patchy[in_shell,4]
    patchy = readdlm(ifilename_r)
    red = patchy[:,5]
    in_shell = (red .< zmax) .& (red .> zmin)
    xyz_3 = patchy[in_shell,1:3]
    w_3 = patchy[in_shell,4]
    for i in 1:length(w_12)
        if w_12[i] != 0
            w_12[i] = 1
        end
    end
    for i in 1:length(w_3)
        if w_3[i] != 0
            w_3[i] = 1
        end
    end
    xyz_12 = transpose(xyz_12)
    w_12 = transpose(w_12)
    xyz_3 = transpose(xyz_3)
    w_3 = transpose(w_3)
    counts = triplet_counts(xyz_12, xyz_3, w_12, w_3, Lsub, 0, rmax, Nbin)
    write_counts(ofilename, counts, 0, 20, 20, sum(w_12))
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end
