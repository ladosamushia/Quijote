using Base.Threads
using Base.Iterators
using DelimitedFiles
using StaticArrays
using Printf

"""
    unique_neighbouring_cubes!(i_n)

Return 14 unique neighbours to a subcube in 3D (including itself).

i_n - must be 3x14.

14 (not 27 because I am only moving in one direction to avoid double counting in a loop.
"""
function unique_neighbouring_cubes!(i_n)
    # Only move forward
    counter = 1
    for di in -1:1, dj in -1:1
        i_n[:,counter] = [di, dj, 1]
        counter += 1
        if dj >= di && dj > -1
            i_n[:,counter] = [di, dj, 0]
            counter += 1
        end
    end
    return nothing
end

"""
    unique_triplet_indeces()

Computes unique triplets with all subcubes adjasent.

This always returns the same 6x71 array.
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

"""
    get_triplet_indeces(r12, r13, r23, dr)

Return three integer indeces randked from least to greatest.

Assumes rmin = 0.
"""
function get_triplet_indeces(r12, r13, r23, dr) 
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
    return imin, imid, imax
end

function double_bin(xyz1, xyz2, w1, w2, dr, rmax, counts)
#    println("started double_bin")
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
#    println("started tri_bin")
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
                imin, imid, imax = get_triplet_indeces(r12, r13, r23, dr)
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
                    i2 = i_n[1,cc] + i1
                    j2 = i_n[2,cc] + j1
                    k2 = i_n[3,cc] + k1
                    # Make sure we are inside the cube
                    if max(i2, j2, k2) <= Nsub && min(i2, j2, k2) >= 1
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
make_cube(xyz, w, Nsub, Lsub)

Take xyz, w vectors and place them into a volume with Nsub subcubes of size Lsub.

xyz must be shifted to zero.
"""
function make_cube(xyz, w, Nsub, Lsub)
    println("started make_cube")
    xyz_cube = Array{Array{SVector{3,Float64}}}(undef, Nsub, Nsub, Nsub)
    w_cube = Array{Array{Float64,1}}(undef, Nsub, Nsub, Nsub)
    for index in eachindex(w)
        i, j, k = ceil.(Int, xyz[:,index]/Lsub)
        if i == 0 
            i = 1
        end
        if j == 0
            j = 1
        end
        if k == 0
            k = 1
        end
        if isassigned(xyz_cube, i, j, k)
            push!(xyz_cube[i, j, k], xyz[:,index])
            push!(w_cube[i, j, k], w[index])
        else
            xyz_cube[i, j, k] = [xyz[:,index],]
            w_cube[i, j, k] = [w[index],]
        end
    end
    return xyz_cube, w_cube
end

"""
write_triplet_counts(ofilename, counts, rmin, rmax, Nbin, Nwgal)

Write triplet counts to output file.
"""
function write_triplet_counts(ofilename, counts, rmin, rmax, Nbin, Nwgal)
    println("started write_triplet_counts")
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

"""
write_pair_counts(ofilename, counts, rmin, rmax, Nbin, Nwgal)

Write pair counts to outpur file.
"""
function write_pair_counts(ofilename, counts, rmin, rmax, Nbin, Nwgal)
    println("started write_pari_counts")
    dr = (rmax - rmin)/Nbin
    ofile = open(ofilename, "w")
    @printf(ofile, "# Weighted number of galaxies: %.1lf\n", Nwgal)
    for i in 1:Nbin
        @printf(ofile, "%.1lf %.1lf\n", dr*(i - 0.5),  counts[1,i])
    end
    close(ofile)
    return nothing
end

"""
read_Patchy(ifilename, zmin, zmax)

Read x, y, z, and weights from a patchy file between redshifts zmin and zmax.
This moves galaxies to xyz>=0.
"""
function read_Patchy(ifilename, zmin, zmax)
    patchy = readdlm(ifilename)
    red = patchy[:,5]
    # Redshift selection
    in_shell = (red .< zmax) .& (red .> zmin)
    xyz = patchy[in_shell,1:3]
    w = patchy[in_shell,4]
    # Only veto weights
    for i in 1:length(w)
        if w[i] != 0
            w[i] = 1
        end
    end
    xyz = transpose(xyz)
    w = transpose(w)
    xyz_min = minimum(xyz, dims=2)
    println("xyz min: ", xyz_min)
    xyz = xyz .- xyz_min
    return xyz, w
end

function dd_count_Patchy(ifile, ofile, zmin, zmax)
    xyz, w = read_Patchy(ifile, zmin, zmax)
    Nw = sum(w)[1]
    xyz_max = maximum(xyz, dims=2)
    L = maximum(xyz_max - xyz_min)
    xyz = xyz .- xyz_min
    Nsub = ceil(Int, L/Lsub)
    dr = rmax/Nbin
    xyz_cube, w_cube = make_cube(xyz, w, Nsub, Lsub)
    pair_hist = zeros(nthreads(), Nbin)
    cube_pairs(xyz_cube, xyz_cube, w_cube, w_cube, Nsub, dr, rmax, pair_hist)
    pair_hist = sum(pair_hist, dims=1)
    write_pair_counts(ofile, pair_hist, 0, rmax, Nbin, Nw)
end

function dr_count_Patchy(dfile, rfile, ofile, zmin, zmax)
    xyz_r, w_r = read_Patchy(rfile, zmin, zmax)
    Nw_r = sum(w_r)[1]
    xyz_max = maximum(xyz_r, dims=2)
    L = maximum(xyz_max - xyz_min)
    xyz_r = xyz_r .- xyz_min
    Nsub = ceil(Int, L/Lsub)
    dr = rmax/Nbin
    xyz_r_cube, w_r_cube = make_cube(xyz_r, w_r, Nsub, Lsub)
    xyz_d, w_d = read_Patchy(dfile, zmin, zmax)
    xyz_d = xyz_d .- xyz_min
    xyz_d_cube, w_d_cube = make_cube(xyz_d, w_d, Nsub, Lsub)
    pair_hist = zeros(nthreads(), Nbin)
    cube_pairs(xyz_r_cube, xyz_d_cube, w_r_cube, w_d_cube, Nsub, dr, rmax, pair_hist)
    pair_hist = sum(pair_hist, dims=1)
    write_pair_counts(ofile, pair_hist, 0, rmax, Nbin, Nw_r)
end

function ddd_count_Patchy(ifilename, ofilename, Lsub, rmax, Nbin, zmin, zmax)
    xyz, w = read_Patchy(ifile, zmin, zmax)
    Nw = sum(w)[1]
    xyz_max = maximum(xyz, dims=2)
    L = maximum(xyz_max - xyz_min)
    xyz = xyz .- xyz_min
    Nsub = ceil(Int, L/Lsub)
    dr = rmax/Nbin
    xyz_cube, w_cube = make_cube(xyz, w, Nsub, Lsub)
    triplet_hist = zeros(nthreads(), Nbin, Nbin, Nbin)
    cube_triplets(xyz_cube, xyz_cube, w_cube, w_cube, Nsub, dr, rmax, triplet_hist) 
    write_counts(ofilename, triplet_hist, 0, rmax, Nbin, Nw)
    return nothing
end

function ddr_count_Patchy(ifilename_d, ifilename_r, ofilename, Lsub, rmax, Nbin, zmin, zmax)
    xyz_r, w_r = read_Patchy(rfile, zmin, zmax)
    Nw_r = sum(w_r)[1]
    xyz_max = maximum(xyz_r, dims=2)
    L = maximum(xyz_max - xyz_min)
    xyz_r = xyz_r .- xyz_min
    Nsub = ceil(Int, L/Lsub)
    dr = rmax/Nbin
    xyz_r_cube, w_r_cube = make_cube(xyz_r, w_r, Nsub, Lsub)
    xyz_d, w_d = read_Patchy(dfile, zmin, zmax)
    xyz_d = xyz_d .- xyz_min
    xyz_d_cube, w_d_cube = make_cube(xyz_d, w_d, Nsub, Lsub)
    xyz_r_cube, w_r_cube = make_cube(xyz_r, w_r, Nsub, Lsub)
    cube_triplets(xyz_d_cube, xyz_r_cube, w_d_cube, w_r_cube, Nsub, dr, rmax, triplet_hist) 
    write_counts(ofilename, triplet_hist, 0, rmax, Nbin, Nw_r)
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end
