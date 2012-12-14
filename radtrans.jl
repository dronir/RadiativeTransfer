#!/usr/bin/env julia

module RadiativeTransfer
export Ray, trace_ray
export point_in_sphere, spherical_start
export point_in_disc, disc_start
export point_in_planar, planar_start

# Radiative transfer in a spherical medium.

# Function for computing unit vector from given vector
unit(v::Vector{Float64}) = v / norm(v)

# Type for ray bookkeeping
type Ray
	origin::Vector{Float64}
	direction::Vector{Float64}
	intensity::Float64
end

# Some functions for checking whether a point is in a medium
point_in_sphere(point::Vector, radius::Real) = norm(point) < radius
point_in_planar(point::Vector) = point[3] < 0
function point_in_disc(point::Vector, radius::Real, depth::Real) 
	return point[3] < depth && norm(point[1:2]) < radius && point[3] > 0
end

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


# Compute the scattered ray given incoming ray, scattering location, asymmetry parameter
# for the HG phase function, and the single-scattering albedo.
function scattered_ray(ray::Ray, location::Vector{Float64}, g::Float64, omega::Float64)
	phi = 2pi * rand()
	theta = random_scattering_angle(g)
	sp = sin(phi)
	cp = cos(phi)
	st = sin(theta)
	ct = cos(theta)
	x = st*sp
	y = st*cp
	z = ct
	
	theta_rotation = acos(ray.direction[3])
	ctr = cos(theta_rotation)
	str = sin(theta_rotation)
	phi_rotation = str==0.0 ? acos(ray.direction[1] / str) : 0.0
	cpr = cos(phi_rotation)
	spr = sin(phi_rotation)
	
	a1 = cpr
	a2 = spr
	b1 = ctr
	b2 = str

	t = (b1*x - b2*z)
	new_x = a1*t - a2*y
	new_y = a2*t + a1*y
	new_z = b2*x + b1*z
	
	ray.origin = location
	ray.direction[1] = new_x
	ray.direction[2] = new_y
	ray.direction[3] = new_z
	ray.intensity *= omega
	return ray
end


# Height of a hemisphere at given radius
sphere_height(R::Float64, r::Float64) = sqrt(R^2 - r^2)

# Function to generate ray starting position in the spherical case
function spherical_start(radius::Real)
	phi = 2pi*rand()
	r = radius*sqrt(rand())
	height = sphere_height(radius, r)
	return Ray([r*cos(phi), r*sin(phi), height], [0.0, 0.0, -1.0], 1.0)
end

# Function to generate ray starting position in disc case
function disc_start(radius::Real, incidence_angle::Real)
	phi = 2pi*rand()
	r = radius*sqrt(rand())
	return Ray([r*cos(phi), r*sin(phi), 0.0], [-sin(incidence_angle), 0.0, -cos(incidence_angle)], 1.0)
end

# Function to generate ray starting position in plane-parallel case
function planar_start(incidence_angle::Real)
	return Ray([0.0, 0.0, 0.0], [-sin(incidence_angle), 0.0, -cos(incidence_angle)], 1.0)
end


# Function to trace one ray into the medium, starting with intensity 1.0
# incoming from the z-direction. Returns the escaping ray.
function trace_ray(starting_ray::Function, inside_medium::Function, max_order::Int64, omega::Float64, g::Float64)
	ray = starting_ray()
	
	# Start the raytracing
	for i = 1:max_order
		tau = random_depth()
		location = ray.origin + tau*ray.direction
		if inside_medium(location)
			ray = scattered_ray(ray, location, g, omega)
		else
			return ray
		end
	end
	return Ray([0.0,0.0,0.0],[0.0,0.0,0.0],0.0)
end

trace_ray(startf::Function, inside::Function, w::Float64, g::Float64) = trace_ray(startf,inside,1000,w,g)
trace_ray(w::Float64, g::Float64) = trace_ray(()->planar_start(0.0), point_in_planar, int(-log(1/w,eps())), w, g)


end # Module