#!/usr/bin/env julia
require("radtrans.jl")
using RadiativeTransfer


function main(N::Integer, omega::Real, phase_params::Vector{Float64})
	detectors = Array(Detector, 3)
	detectors[1] = Detector([0.0, 0.0, 1.0])
    detectors[2] = Detector(unit([1.0, 0.0, 1.0]))
    detectors[3] = Detector(unit([-1.0, 0.0, 1.0]))
    
    
	startf() = planar_start(.785398163)
	medium = PlanarMedium()
	
	println("\nTracing $N rays...")
	tic()
	for i = 1:N
		detectors = trace_ray(detectors, startf, medium, 10, omega, phase_params)
	end
	toc()
	
	println(detectors)
end

