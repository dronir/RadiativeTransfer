#!/usr/bin/env julia
require("radtrans.jl")
using RadiativeTransfer

function main(N::Integer)
	scattered_rays = Array(Ray, N)

	println("\nTracing $N rays...")
	tic()
	for i = 1:N
		scattered_rays[i] = trace_ray(2.0, 0.9, 0.6)
	end
	toc()
	
	scattered_intensity = zeros(31)
	scattering_angles = zeros(N)
	I_scattered = 0.0
	for i = 1:N
		ray = scattered_rays[i]
		I_scattered += ray.intensity
		theta = acos(dot(ray.direction, [0.0, 0.0, 1.0]))
		scattering_angles[i] = theta
		exit_bin = int((theta / pi) * 30) + 1
		scattered_intensity[exit_bin] += ray.intensity
	end
	

	@printf("  Scattered intensity: %.4f\n", I_scattered/N)
	@printf("   Absorbed intensity: %.4f\n", (1 - I_scattered/N))
	@printf("  Mean scattering angle: %.2f°\n", mean(scattering_angles)*180/pi)
	csvwrite("angles.txt", scattering_angles)
	println("Wrote scattering angles to angles.txt.")
	csvwrite("intensity.txt", scattered_intensity)
	println("Wrote scattering intensity to intensity.txt.")
	println("")
end

main(10000)
