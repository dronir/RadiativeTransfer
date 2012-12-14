#!/usr/bin/env julia
require("radtrans.jl")
using RadiativeTransfer

function main(N::Integer, tau::Float64, omega::Float64, g::Float64)
	scattered_rays = Array(Ray, N)

	omega = 0.9

	println("\nTracing $N rays...")
	tic()
	for i = 1:N
		scattered_rays[i] = trace_ray(tau, omega, g)
	end
	toc()
	
	scattered_intensity = zeros(31)
	scattering_angles = zeros(N)
	scattering_orders = zeros(N)
	I_scattered = 0.0
	for i = 1:N
		ray = scattered_rays[i]
		I_scattered += ray.intensity
		theta = acos(dot(ray.direction, [0.0, 0.0, 1.0]))
		scattering_angles[i] = theta
		exit_bin = int(floor((theta / pi) * 30)) + 1
		exit_bin = min(exit_bin, 30)
		low = ((exit_bin-1)*pi/30)
		hi = ((exit_bin)*pi/30)
		A = 2pi*(cos(low)-cos(hi))
		scattered_intensity[exit_bin] += ray.intensity / A
		scattering_orders[i] = ceil(-1/log(scattered_rays[i].intensity, 2))
	end
	

	@printf("  Scattered intensity: %.4f\n", I_scattered/N)
	@printf("   Absorbed intensity: %.4f\n", (1 - I_scattered/N))
	@printf("  Mean scattering angle: %.2f°\n", mean(scattering_angles)*180/pi)
	csvwrite("angles.txt", scattering_angles)
	println("Wrote scattering angles to angles.txt.")
	csvwrite("intensity.txt", scattered_intensity)
	csvwrite("orders.txt", scattering_orders)
	println("Wrote scattering intensity to intensity.txt.")
	println("")
end

main(N::Integer) = main(N, 2.0, 0.9, 0.5)

#main(10000)
