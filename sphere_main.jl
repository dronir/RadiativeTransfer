#!/usr/bin/env julia
require("radtrans.jl")
using RadiativeTransfer




function trace_geometry(N::Integer, alpha::Real, omega::Real, phase_params::Vector{Float64})
    radius = 20.0
    max_order = 20
	detectors = Array(Detector, 1)
	detectors[1] = Detector([sin(alpha), 0.0, cos(alpha)], max_order)
    
	startf() = spherical_start(radius)
	medium = SphericalMedium(radius)
	
	for i = 1:N
		trace_ray(detectors, startf, medium, max_order, omega, phase_params)
	end
	
	return detectors[1].intensity
end


function main(N::Integer, M::Integer, omega::Real, phase_params::Vector{Float64})
    alphas = linspace(0.0+eps(), pi-eps(), M)
    traced_first_order = zeros(M)
    traced_total = zeros(M)
    analytical_LS = zeros(M)
    tic()
    for i in 1:M
        alpha = alphas[i]
        orders = trace_geometry(N, alpha, omega, phase_params)
        traced_first_order[i] = orders[1]
        traced_total[i] = sum(orders)
        analytical_LS[i] = pi/32 * omega * double_HG(pi-alpha, phase_params) * (1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4)))
    end
    toc()
    traced_first_order /= traced_first_order[1]
    traced_total /= traced_total[1]
    analytical_LS /= analytical_LS[1]
    results = hcat(alphas, traced_first_order, traced_total, analytical_LS)
    writecsv("output.txt", results)
end
