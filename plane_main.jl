#!/usr/bin/env julia
require("radtrans.jl")
using RadiativeTransfer


function main(N::Integer, theta_i::Real, omega::Real, phase_params::Vector{Float64})
    max_order = 10
    detectors = Array(Detector, 3)
    detectors[1] = Detector([0.0, 0.0, 1.0], max_order)
    detectors[2] = Detector(unit([1.0, 0.0, 1.0]), max_order)
    detectors[3] = Detector(unit([-1.0, 0.0, 1.0]), max_order)
    
    
    startf() = planar_start(theta_i)
    medium = PlanarMedium()
    
    println("\nTracing $N rays...")
    tic()
    for i = 1:N
        trace_ray(detectors, startf, medium, max_order, omega, phase_params)
    end
    toc()
    
    for d in detectors
        println(d.intensity)
        println(sum(d.intensity))
    end
end

