function read_binary_array(N, file)
    X = Vector{UInt8}(undef, 496)
    readbytes!(file, X, N*4)
    X = reinterpret(Float32, N*4)
end
# snapfile = "C:/Users/Lado/Documents/Globus/snap_002.0"
# Ngal = 16596561
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

#Vin_bin = zeros(Int32, 100, 100)

function get_vinz(Ngal, X, Y, Z, VZ, distmax, VinZmin, VinZmax)
    Vin_bin = zeros(Int32, nthreads(), 100, 100)
    @threads for i in 1:Ngal
        for j in i:Ngal
            dist = trunc(Int32, sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2 + (Z[i] - Z[j])^2)) + 1
            if dist > distmax
                continue
            end
            VinZ = trunc(Int32, (VZ[i] - VZ[j])*sign(Z[i] - Z[j])) + 1
            if VinZ > VinZmax || VinZ < VinZmin
                continue
            end
            Vin_bin[threadid(), dist, VinZ] += 1
        end
    end
    return Vin_bin
end
