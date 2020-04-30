"""
    unique_neighbouring_cubes(i_n)

Return 14 unique neighbours to a subcube in 3D (including itself).

i_n is 3x14.

14 (not 27 because I am only moving in one direction to avoid double counting in a loop.
"""
function unique_neighbouring_cubes()
    i_n = zeros(Int, 3, 14)
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
    return i_n
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

