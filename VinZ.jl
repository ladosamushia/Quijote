using DelimitedFiles
using Base.Threads

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

function get_vinz(Ngal, X, Y, Z, VZ, distmax, VinZmin, VinZmax, L)
    Vin_bin = zeros(Int32, nthreads(), 100, 200)
    dx::Float32 = 0.0
    dy::Float32 = 0.0
    dz::Float32 = 0.0
    @threads for i in 1:Ngal
        for j in i:Ngal
            dx = abs(X[i] - X[j])
            if dx > L/2
                dx -= L/2
            end
            if dx > distmax
                continue
            end
            dy = abs(Y[i] - Y[j])
            if dy > L/2
                dy -= L/2
            end
            if dy > distmax
                continue
            end
            dz = abs(Z[i] - Z[j])
            if dy > L/2
                dz -= L/2
            end
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
            Vin_bin[threadid(), dist_bin, VinZ_bin] += 1
        end
    end
    return Vin_bin
end

function main()
    snapfile = "/home/lado/snap_002.0"
    Ngal = 16596561
    X, Y, Z, VX, VY, VZ = load_snap(snapfile, Ngal)
    X = X[1:100:end]
    Y = Y[1:100:end]
    Z = Z[1:100:end]
    VX = VX[1:100:end]/100
    Vin_bin = get_vinz(size(X)[1], Z, Y, X, VX, 100.0, -25.0, 25.0, 1000.0)
    Vin_bin = sum(Vin_bin, dims=1)
    Vin_bin = reshape(Vin_bin, (100, 200))
    writedlm("VinZ.csv", Vin_bin)
end
