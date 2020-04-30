using Base.Threads
using StaticArrays
include("geometry.jl")

"""
    loop_over_pairs(xyz1, xyz2, w1, w2, dr, rmax, histogram, f_bin)

Go over pairs of arrays and compute counts based on f_bin.
"""
function loop_over_pairs(xyz1, xyz2, w1, w2, dr, rmax, histogram, f_bin)
    for i1 in 1:length(xyz1)
    # if it is a self count do not double count
        if xyz1 == xyz2
            i2min = i1
        else
            i2min = 1
        end
        for i2 in i2min:length(xyz2)
            r12 = sqrt(sum((xyz1[i1] - xyz2[i2]).^2))
            # Enforce maximum distance
            if r12 >= rmax || r12 <= 1e-3
                continue
            end
            f_bin(xyz1[i1], xyz2[i2], w1[i1], w2[i2], dr, histogram)
        end
    end
end

"""
    cube_pairs(xyz_cube_1, xyz_cube_2, w_cube_1, w_cube_2, dr, rmax, histogram, f_bin)

Go through all pairs of cubes, and all pairs of particles in those cubes and compute histogram of counts defined in f_bin)
"""
function cube_pairs(xyz_cube_1, xyz_cube_2, w_cube_1, w_cube_2, dr, rmax, histogram, f_bin)
    println("start cube_pairs")
    i_n = unique_neighbouring_cubes()
    Nsub = size(xyz_cube_1)[1]
    # All subcubes
    @threads for i1 in 1:Nsub
        for j1 in 1:Nsub, k1 in 1:Nsub
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
                        loop_over_pairs(xyz1, xyz2, w1, w2, dr, rmax, histogram, f_bin)
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
function cube_triplets(xyz_cube_12, xyz_cube_3, w_cube_12, w_cube_3, Nsub, dr, rmax, histogram, f_bin)
    println("start cube_triplets")
    j_n = unique_triplet_indeces()
    # All subcubes
    @threads for i1 in 1:Nsub
        for j1 in 1:Nsub, k1 in 1:Nsub
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
                        f_bin(xyz1, xyz2, xyz3, w1, w2, w3, dr, rmax, histogram)
                    end
                end
            end
        end
    end
end

"""
make_cube(xyz, w, Lsub)

Take xyz, w vectors and place them into a volume with subcubes of size Lsub.

xyz must be shifted to zero.
"""
function make_cube(xyz, w, Lsub)
    println("started make_cube")
    xyz_max = maximum(xyz)
    Nsub = ceil.(Int, xyz_max/Lsub)
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
