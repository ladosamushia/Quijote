using DelimitedFiles
using Base.Threads
using StaticArrays
using CUDAdrv, CUDAnative, CuArrays

function read_binary_array(N, file)
    X = Vector{UInt8}(undef, N*4)
    readbytes!(file, X, N*4)
    X = reinterpret(Float32, X)
end

function load_snap(snapfile, Ngal)
    snapfile = open(snapfile, "r")
    # Read header
    temp_array = Vector{UInt8}(undef, 496)
    readbytes!(snapfile, temp_array, 496)
    # Read X, Y, Z
    X = read_binary_array(Ngal, snapfile)
    Y = read_binary_array(Ngal, snapfile)
    Z = read_binary_array(Ngal, snapfile)
    VX = read_binary_array(Ngal, snapfile)
    VY = read_binary_array(Ngal, snapfile)
    VZ = read_binary_array(Ngal, snapfile)

    println(minimum(X), " ", maximum(X))
    println(minimum(Y), " ", maximum(Y))
    println(minimum(Z), " ", maximum(Z))
    println(minimum(VX), " ", maximum(VX))
    println(minimum(VY), " ", maximum(VY))
    println(minimum(VZ), " ", maximum(VZ))

    X /= 1000
    Y /= 1000
    Z /= 1000

    return (X, Y, Z, VX, VY, VZ)
end

function shortest_distance(a, b, L)
    dist = abs(a - b)
    if dist > L/2
        dist == L/2
    end
    return dist
end

function get_vinz(Ngal, X, Y, Z, VZ, distmax, VinZmin, VinZmax, L, Vin_hist, dist_hist)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i in index:stride:Ngal
        for j in i:Ngal
            dx = shortest_distance(X[i], X[j], L)
            if dx > distmax
                continue
            end
            dy = shortest_distance(Y[i], Y[j], L)
            if dy > distmax
                continue
            end
            dz = shortest_distance(Z[i], Z[j], L)
            if dz > distmax
                continue
            end
            dist = sqrt(dx^2 + dy^2 + dz^2)
            if dist >= distmax
                continue
            end
            dist_bin = floor(Int32, dist) + 1
            VinZ = (VZ[j] - VZ[i])*sign(Z[i] - Z[j])
            if VinZ >= VinZmax || VinZ <= VinZmin
                continue
            end
            VinZ_bin = floor(Int32, VinZ*4.0) + 101
            dist_hist[threadIdx().x, dist_bin] += 1
            Vin_hist[threadIdx().x, dist_bin, VinZ_bin] += 1
        end
    end
    return 
end

function main()
    snapfile = "/home/lado/snap_002.0"
    Ngal = 16596561
    Xall, Yall, Zall, VXall, VYall, VZall = load_snap(snapfile, Ngal)
    X = CuArray(Xall[1:10:end])
    Y = CuArray(Yall[1:10:end])
    Z = CuArray(Zall[1:10:end])
    println(typeof(X))
    println(sizeof(X))
    VX = CuArray(VXall[1:10:end]/100)
    Vin_hist = zeros(Int32, 512, 100, 200)
    dist_hist = zeros(Int32, 512, 200)
    Vin_hist = CuArray(Vin_hist)
    dist_hist = CuArray(dist_hist)
    @cuda threads=512 blocks=256 get_vinz(size(X)[1], Z, Y, X, VX, 100.0, -25.0, 25.0, 1000.0, Vin_hist, dist_hist)
    Vin_hist = Array(Vin_hist)
    dist_hist = Array(dist_hist)
    Vin_hist = sum(Vin_hist, dims=1)
    dist_hist = sum(dist_hist, dims=1)
    Vin_hist = reshape(Vin_hist, (100, 200))
    writedlm("VinZ.csv", Vin_hist)
    writedlm("dist.csv", dist_hist)
end

main()
