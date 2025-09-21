using HDF5
using LinearAlgebra
using Makie 
using GLMakie


global const fid::HDF5.File          = h5open("/home/nabil/work/projects/hrfile_hdf5_converter/cmake-build-debug/wann.hdf5", "r") 
global const nw::Int64               = read(fid["nw"])
global const nr::Int64               = read(fid["nr"]);

global const h::Array{ComplexF64, 1} = read(fid["reH"]) + im * read(fid["imH"])
global const r::Matrix{ComplexF64}   = reshape(read(fid["rvecs"]), (3, nr))';

close(fid)


# Lattice Constants 
global const A₁::Vector{Float64} = [3.9413692012334867, -0.0000000000000000, 0.0000000000000000] 
global const A₂::Vector{Float64} = [0.0000000000000000,  3.9413692012334867, 0.0000000000000000]
global const A₃::Vector{Float64} = [0.0000000000000000, 0.0000000000000000,  3.9413692012334867]
global const A::Matrix{Float64}  = hcat(A₁, A₂, A₃)'  
global const G::Matrix{Float64}  = 2π * inv(A)'
global const G₁::Vector{Float64} = G[1, :]
global const G₂::Vector{Float64} = G[2, :]
global const G₃::Vector{Float64} = G[3, :];


global const nk::Int64 = 100; # Number of k-points along each direction

# K-Path 
global const nodes::Vector{Vector{Float64}} = [[0.0, 0.0, 0.0],
         [0.5, 0.5, 0.0],
         [0.5, 0.0, 0.0],
         [0.0, 0.0, 0.0]];

function kpath_connect_nodes_red(nodes, nk)
    ks_red = Vector{Vector{Float64}}()
    for ik ∈ 1 : length(nodes) - 1
        start_node = nodes[ik]
        end_node   = nodes[ik + 1]
        pts        = LinRange(start_node, end_node, nk)
        if ik != length(nodes) - 1
            pts = pts[1 : end - 1]
        end
        append!(ks_red, pts)
    end
    reduce(hcat, ks_red)'
end


# Calculate the k-path in cartesian coordinates 
function red_to_cart(ks_red, G₁, G₂, G₃)
    ks_cart::Matrix{Float64} = zeros(Float64, size(ks_red, 1), 3)
    for ik ∈ axes(ks_red, 1)
        k_red  = ks_red[ik, :]
        ks_cart[ik, :] = k_red[1] * G₁ + k_red[2] * G₂ + k_red[3] * G₃
    end
    ks_cart
end


function calc_xs(ks_cart)
    xs = Vector{Float64}(undef, size(ks_cart, 1))
    xs[1] = 0.0
    for ik ∈ 2 : size(ks_cart, 1)
        dk = ks_cart[ik, :] - ks_cart[ik - 1, :]
        xs[ik] = xs[ik - 1] + norm(dk)
    end
    xs
end

global const ks_red::Matrix{Float64}  = kpath_connect_nodes_red(nodes, nk);
global const ks_cart::Matrix{Float64} = red_to_cart(ks_red, G₁, G₂, G₃);
global const xs::Vector{Float64}      = calc_xs(ks_cart);


# Compute e^ikR Matrix 
function calculate_eikR_mat(ks_red::Matrix{Float64}, Rs::Matrix{ComplexF64}, nk, nr)
    eikR_mat::Matrix{ComplexF64} = zeros(ComplexF64, nk, nr)
    @inbounds Threads.@threads for ik ∈ 1 : nk
        @inbounds for iR ∈ 1 : nr
            eikR_mat[ik, iR] = cispi(2 * dot(ks_red[ik, :], Rs[iR, :]))
        end
    end
    eikR_mat
end

@inline function ek(h, eikR, nw, nr)
    Res::Matrix{Float64}     = zeros(Float64, nw, size(eikR, 1))
    for ik ∈ 1 : size(eikR, 1)
        Rsum::Matrix{ComplexF64} = zeros(ComplexF64, nw, nw)
        for α ∈ 0 : nw - 1 
            for β ∈ 0 : nw - 1 
                for iR ∈ 0 : size(eikR, 2) - 1
                    Rsum[α + 1, β + 1] += h[iR + nr * nw * α + nr * β + 1] * eikR[ik, iR + 1]
                end # β
            end   # α
        end # iR
        Res[:, ik] = eigvals(Hermitian(Rsum))
    end # ik 
        
    return Res
end

global const eikR_mat::Matrix{ComplexF64} = calculate_eikR_mat(ks_red, r, size(ks_red, 1), nr);
global const eks::Matrix{Float64}         = ek(h, eikR_mat, nw, nr);


let 
    f = Figure(size = (1000, 1000))
    ax = Axis(f[1, 1])
    for i ∈ 1 : nw
        lines!(ax, xs, R[i, :], color = :blue)
    end
    save("bands.png", f)
end


