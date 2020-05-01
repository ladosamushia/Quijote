using Base.Threads
using StaticArrays
include("geometry.jl")
include("io.jl")
include("cubing.jl")

"""
    distance_bin(xyz1, xyz2, w1, w2, dr, histogram)

Weighted histogram of pairs distances between two points.
"""
function distance_bin(xyz1, xyz2, w1, w2, dr, histogram)
    r12 = sqrt(sum((xyz1 - xyz2).^2))
    index = ceil(Int, r12/dr)
    histogram[threadid(),index] += w1*w2
    return nothing
end

"""
    wp_bin(xyz1, xyz2, w1, w2, dr, histogram)

Weighted histogram of angular distances.
"""
function wp_bin(xyz1, xyz2, w1, w2, dr, histogram)
    r12 = sqrt(sum((xyz1[1:2] - xyz2[1:2]).^2))
    index = ceil(Int, r12/dr)
    histogram[threadid(),index] += w1*w2
    return nothing
end

"""
    wp_survey_bin(xyz1, xyz2, w1, w2, dr, histogram)

Weighted histogram of angular distances but accounting for varying LOS.
"""
function wp_survey_bin(xyz1, xyz2, w1, w2, dr, histogram)
    r1 = sqrt(sum(xyz1.^2))
    r2 = sqrt(sum(xyz2.^2))
    if r2 > r1
        r12 = sqrt(sum((xyz2*r1/r2 - xyz1).^2))
    else
        r12 = sqrt(sum((xyz1*r2/r1 - xyz2).^2))
    end
    index = ceil(Int, r12/dr)
    histogram[threadid(),index] += w1*w2
    return nothing
end

"""
    vin_bin(xyz1, xyz2, v1, v2, dr, histogram)

Weighted histogram of infall velocities.
"""
function vin_bin(xyz1, xyz2, v1, v2, dr, histogram)
    r12 = sqrt(sum((xyz1 - xyz2).^2))
    if xyz1[3] < xyz2[3]
        v12 = v1 - v2
    else
        v12 = v2 - v1
    end
    index_v = ceil(Int, (v12 + 51)/100)
    if index_v < 1 || v12 > 100
        return nothing
    end
    index_r = ceil(Int, r12/dr)
    histogram[threadid(),index_r,index_v] += 1
    return nothing
end

"""
    tri_bin(xyz1, xyz2, xyz3, w1, w2, w3, dr, rmax, histogram)

bin distances between particles in three arrays and incriment histogram in counts.
"""
function tri_bin(xyz1, xyz2, xyz3, w1, w2, w3, dr, rmax, histogram)
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
                histogram[threadid(),imin,imid,imax] += w1[i1]*w2[i2]*w3[i3]
            end
        end
    end
    return nothing
end
