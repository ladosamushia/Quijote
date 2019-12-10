using Base.Threads
using Base.Iterators
using DelimitedFiles

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
    for i = 1:27
        for j = 1:3
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
        # if xyzw1 == xyzw2
            # i2min = i1
        # else
            # i2min = 1
        # end
        for i2 in 1:length(xyzw2)
            r12 = sqrt(sum((xyzw1[i1][1:3] - xyzw2[i2][1:3]).^2))
            if r12 >= rmax || r12 == 0
                continue
            end
            # if xyzw2 == xyzw3
                # i3min = i2
            # else
                # i3min = 1
            # end
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

    # All cubes
    index = []
    for i in 1:Nsub, j in 1:Nsub, k in 1:Nsub
        push!(index, [i, j, k])
    end
    Ncubes = length(index)
    # All unique cube triplets
    tri_index = []
    for i in 1:Ncubes, j in i:Ncubes, k in j:Ncubes
        # Only keep the neighbours
        if maximum(abs.(index[i] - index[j])) <= 1 && maximum(abs.(index[k] - index[j])) <=1 && maximum(abs.(index[i] - index[k])) <=1
            push!(tri_index, [index[i], index[j], index[k]])
        end
    end
    
    @threads for index_tri in tri_index
        ijk1, ijk2, ijk3 = index_tri
        xyzw1 = xyzw_cube[ijk1[1],ijk1[2],ijk1[3]]
        xyzw2 = xyzw_cube[ijk2[1],ijk2[2],ijk2[3]]
        xyzw3 = xyzw_cube[ijk3[1],ijk3[2],ijk3[3]]
        tri_bin(xyzw1, xyzw2, xyzw3, dr, rmax, counts)
    end 

end

function triplet_counts(Ngal, Lsurvey, xyzw, Lsub, rmin, rmax, Nbin)

    println("Ngal ", Ngal)
    Nsub = ceil(Int, Lsurvey/Lsub)
    println("Nsub ", Nsub)
    dr = (rmax - rmin)/Nbin
    println(dr)
    println(nthreads())
    counts = zeros(Float32, nthreads(), Nbin, Nbin, Nbin)
    println(typeof(counts))
    # Fill in the sub-cubes
    xyzw_cube = Array{Array{Array{Float64,1}}}(undef, Nsub, Nsub, Nsub)
    min_xyz = [minimum(xyzw[1,:]), minimum(xyzw[2,:]), minimum(xyzw[3,:])]
    i_xyz = ceil.(Int, (xyzw[1:3,:] .- min_xyz)/Lsub)
    for i in 1:Ngal
        if i_xyz[1,i] == 0
            i_xyz[1,i] = 1
        end
        if i_xyz[2,i] == 0
            i_xyz[2,i] = 1
        end
        if i_xyz[3,i] == 0
            i_xyz[3,i] = 1
        end
    end
    for i = 1:Nsub, j = 1:Nsub, k = 1:Nsub
        xyzw_cube[i, j, k] = [zeros(4),]
    end
    for i in 1:Ngal
        push!(xyzw_cube[i_xyz[1,i],i_xyz[2,i],i_xyz[3,i]], xyzw[:,i])
    end
    # println(xyzw_cube)
    # Count triplets
    println("cube_triplets")
    cube_triplets(xyzw_cube, Nsub, dr, rmax, counts)
    counts = sum(counts, dims=1)
    println("sum", sum(counts))
    println(size(counts))
    return counts

end

function test_triplet_counts()
    xyzw = rand(4, 1000000)
    xyzw[1:3,:] *= 1000
    triplet_counts(1000000, 1000, xyzw, 20, 0, 20, 20)
    return nothing
end

function test_on_Patchy()
    patchy = readdlm("Patchy-Mocks-DR12NGC-COMPSAM_V6C_0001_xyzw.dat")
    xyzw = patchy[:,1:4]
    xyzw[:,4] .= 1
    xyzw = transpose(xyzw)
    Ngal = size(xyzw)[2]
    Lsurvey = 3300
    triplet_counts(Ngal, Lsurvey, xyzw, 400, 0, 20, 20)
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end
