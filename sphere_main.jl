#!/usr/bin/env julia

using RadiativeTransfer

# Trace N rays to detectors at phase angles 'alphas' using given omega and 2HG parameters
function trace(N::Integer, alphas::Array{Float64}, omega::Real, phase_params::Vector{Float64})
    radius = 20.0
    max_order = 20
    M = length(alphas)
    
    # Put a detector at each phase angle
    detectors = Array(Detector, M)
    for i = 1:M
        alpha = alphas[i]
        detectors[i] = Detector([sin(alpha), 0.0, cos(alpha)], max_order)
    end
    
    # Set up medium
    startf() = spherical_start(radius)
    medium = SphericalMedium(radius)
    
    # Trace N rays
    for i = 1:N
        trace_ray(detectors, startf, medium, max_order, omega, phase_params)
    end
    return detectors
end

# Simulate Lommel-Seeliger with N rays, M steps of alpha with given omega and 2HG parameters
function main(N::Integer, M::Integer, omega::Real, phase_params::Vector{Float64})
    
    alphas = linspace(0.0+eps(), pi-eps(), M)
    traced_first_order = zeros(M)
    traced_total = zeros(M)
    analytical_LS = zeros(M)
    
    # Do the raytracing
    tic()
    detectors = trace(N, alphas, omega, phase_params)
    toc()
    
    # Extract the intensities from the detectors and compute analytical values
    for i in 1:M
        alpha = alphas[i]
        orders = detectors[i].intensity
        traced_first_order[i] = orders[1]
        traced_total[i] = sum(orders)
        analytical_LS[i] = pi/32 * omega * double_HG(pi-alpha, phase_params) * (1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4)))
    end
    
    # Normalize values to zero phase angle
    traced_total /= traced_first_order[1]
    traced_first_order /= traced_first_order[1]
    analytical_LS /= analytical_LS[1]
    
    # Save results to file
    results = hcat(alphas, traced_first_order, traced_total, analytical_LS)
    writecsv("output.txt", results)
end
