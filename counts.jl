include("triplet_count.jl")
include("cubing.jl")

"""
    dd_count_Patchy(ifile, ofile, zmin, zmax, rmax, Nbin)
    
compute and write to a file dd counts from Patchy.
"""
function dd_count_Patchy(ifile, ofile, zmin, zmax, rmax, Nbin)
    xyz, w = read_Patchy(ifile, zmin, zmax)
    xyz_cube, w_cube = make_cube(xyz, w, rmax)
    pair_hist = zeros(nthreads(), Nbin)
    dr = rmax/Nbin
    cube_pairs(xyz_cube, xyz_cube, w_cube, w_cube, dr, rmax, pair_hist, wp_survey_bin)
    pair_hist = sum(pair_hist, dims=1)
    write_pair_counts(ofile, pair_hist, 0, rmax, Nbin, sum(w))
end

function dr_count_Patchy(dfile, rfile, ofile, zmin, zmax, rmax, Nbin)
    xyz_d, w_d = read_Patchy(dfile, zmin, zmax)
    xyz_d_cube, w_d_cube = make_cube(xyz_d, w_d, rmax)
    xyz_r, w_r = read_Patchy(rfile, zmin, zmax)
    xyz_r_cube, w_r_cube = make_cube(xyz_r, w_r, rmax)
    pair_hist = zeros(nthreads(), Nbin)
    dr = rmax/Nbin
    cube_pairs(xyz_d_cube, xyz_r_cube, w_d_cube, w_r_cube, dr, rmax, pair_hist, wp_survey_bin)
    pair_hist = sum(pair_hist, dims=1)
    write_pair_counts(ofile, pair_hist, 0, rmax, Nbin, sum(w))
end

function ddd_count_Patchy(dfile, rfile, ofile, Lsub, rmax, Nbin, zmin, zmax)
    xyz_cube, w_cube, Nw = read_Patchy(dfile, zmin, zmax, Lsub)
    triplet_hist = zeros(nthreads(), Nbin, Nbin, Nbin)
    cube_triplets(xyz_cube, xyz_cube, w_cube, w_cube, Nsub, dr, rmax, triplet_hist, tri_bin)
    triplet_hist = sum(triplet_hist, dims=1)
    write_counts(ofilename, triplet_hist, 0, rmax, Nbin, Nw)
    return nothing
end

function ddr_count_Patchy(dfile, rfile, ofile, Lsub, rmax, Nbin, zmin, zmax)
    xyz_r_cube, w_r_cube, Nr = read_Patchy(rfile, zmin, zmax, Lsub)
    xyz_d_cube, w_d_cube, Nd = read_Patchy(dfile, zmin, zmax, Lsub)
    triplet_hist = zeros(nthreads(), Nbin, Nbin, Nbin)
    cube_triplets(xyz_d_cube, xyz_r_cube, w_d_cube, w_r_cube, Nsub, dr, rmax, triplet_hist, tri_bin)
    triplet_hist = sum(triplet_hist, dims=1)
    write_counts(ofile, triplet_hist, 0, rmax, Nbin, Nd)
    return nothing
end

function three_pcf(DDD, DDR, DRR, RRR, Ngal, Nran)
    alpha = RRR/DDD
    tpcf = (alpha*alpha*alpha*DDD - 3*alpha*alpha*DDR + 3*alpha*DRR - RRR)./RRR
    return tpcf
end
