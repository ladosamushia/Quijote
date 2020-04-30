using Printf
using DelimitedFiles

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
function read_Patchy(ifile, zmin, zmax)
    patchy = readdlm(ifile)
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
    xyz_max = maximum(xyz, dims=2)
    xyz_min = minimum(xyz, dims=2)
    L = maximum(xyz_max - xyz_min)
    println("xyz min: ", xyz_min)
    xyz = xyz .- xyz_min
    println("xyz new min: ", minimum(xyz, dims=2))
    return xyz, w
end
