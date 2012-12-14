#!/usr/bin/env julia

module RadiativeTransfer
export Ray, trace_ray

# Radiative transfer in a spherical medium.

# Function for computing unit vector from given vector
unit(v::Vector{Float64}) = v / norm(v)

# Type for ray bookkeeping
type Ray
	origin::Vector{Float64}
	direction::Vector{Float64}
	intensity::Float64
end

# Some functions related to spheres
point_in_sphere(point::Vector, radius::Real) = norm(point) < radius
sphere_height(R::Float64, r::Float64) = sqrt(R^2 - r^2)

# Generate random free path
random_depth() = -log(rand())

# Generate a random direction weighed with a Henyey-Greenstein
# phase function.
function random_scattering_angle(g::Float64)
	s = 2*rand() - 1
	if g != 0.0
		mu = (1 + g^2 - ((1-g^2)/(1+g*s))^2) / (2g)
	else
		mu = s
	end
	return acos(mu)
end

# Two vector rotations
function rotate_Z(v::Vector{Float64}, theta::Float64)
	ct = cos(theta)
	st = sin(theta)
	M = [[ct, st, 0] [-st, ct, 0] [0, 0, 1]]
	return M * v
end

function rotate_Y(v::Vector{Float64}, theta::Float64)
	ct = cos(theta)
	st = sin(theta)
	M = [[ct, 0, st] [0, 1, 0] [-st, 0, ct]]
	return M * v
end

# Compute the scattered ray given incoming ray, scattering location, asymmetry parameter
# for the HG phase function, and the single-scattering albedo.
function scattered_ray(ray::Ray, location::Vector{Float64}, g::Float64, omega::Float64)
	phi = 2pi * rand()
	theta = random_scattering_angle(g)
	sp = sin(phi)
	cp = cos(phi)
	st = sin(theta)
	ct = cos(theta)
	random_direction = [st*sp, st*cp, ct]
	
	theta_rotation = acos(dot(ray.direction, [0.0, 0.0, 1.0]))
	phi_rotation = acos(ray.direction[1] / (norm(ray.direction[1:2])))
	phi_rotation = isnan(phi_rotation) ? 0.0 : phi_rotation
	
	scattering_direction = rotate_Z(rotate_Y(random_direction, -theta_rotation), -phi_rotation)
	scattered_intensity = ray.intensity * omega
	return Ray(location, scattering_direction, scattered_intensity)
end

# Function to trace one ray into the medium, starting with intensity 1.0
# incoming from the z-direction. Returns the escaping ray.
function trace_ray(tau_R::Float64, omega::Float64, g::Float64)
	phi = 2pi*rand()
	r = tau_R*sqrt(rand())
	height = sphere_height(tau_R, r)
	
	ray = Ray([r*cos(phi), r*sin(phi), height], [0.0, 0.0, -1.0], 1.0)
	
	# Start the raytracing
	while true
		tau = random_depth()
		location = ray.origin + tau*ray.direction
		if point_in_sphere(location, tau_R)
			ray = scattered_ray(ray, location, g, omega)
		else
			return ray
		end
	end
end

end # Module