#!/usr/bin/env julia

# Radiative transfer in a spherical medium.
module RadiativeTransfer


export Ray, trace_ray
export spherical_start, disc_start, planar_start
export Medium, SphericalMedium, PlanarMedium, CylindricalMedium
export Detector
export unit
export double_HG

# Function for computing unit vector from given vector
unit(v::Vector{Float64}) = v / norm(v)

# RAY
# Type for ray bookkeeping
type Ray
    origin::Vector{Float64}
    direction::Vector{Float64}
    intensity::Float64
end

## DETECTOR
# Detector at infinity (to count intensity emitted to that direction)
type Detector
    direction::Vector{Float64}
    intensity::Vector{Float64}
end
# Constructor with intensity count zero
Detector(direction::Vector{Float64}, order::Integer) = Detector(direction, zeros(order))


# MEDIUM
# Define some types for the medium
abstract Medium

# z=0 plane
immutable PlanarMedium <: Medium
end

# Sphere with given radius centered at origin
immutable SphericalMedium <: Medium
    radius::Float64
end

# Cylinder with given radius and half-thickness, centered at origin
immutable CylindricalMedium <: Medium
    radius::Float64
    half_thickness::Float64
end

# Function to check if a given point is inside a medium
point_in_medium(M::SphericalMedium, point::Vector{Float64}) = norm(point) < M.radius
point_in_medium(M::PlanarMedium, point::Vector{Float64}) = point[3] < 0
function point_in_medium(M::CylindricalMedium, point::Vector{Float64})
    return point[3] < M.half_thickness && point[3] > -M.half_thickness && norm(point[1:2]) < M.radius
end

# Get distance in semi-infinite planar medium in given direction from given point
function trace_to_direction(M::PlanarMedium, point::Vector{Float64}, direction::Vector{Float64})
    if direction[3] <= 0
        return Inf
    end
    return -point[3] / direction[3]
end

# Get distance in spherical medium in given direction from given point
function trace_to_direction(M::SphericalMedium, point::Vector{Float64}, direction::Vector{Float64})
    b = 2*dot(direction, point)
    c = dot(point,point) - M.radius^2
    q = -0.5 * (b - sign(b) * sqrt(b^2 - 4*c))
    return sign(q) > 0 ? q : c/q
end

# Get distance in cylindrical medium in given direction from given point
function trace_to_direction(M::CylindricalMedium, point::Vector{Float64}, direction::Vector{Float64})
    dz = direction[3]
    if dz > 0
        t_plane = (point[3] - M.half_thickness) / dz
    elseif dz < 0
        t_plane = -(point[3] + M.half_thickness) / dz
    else
        t_plane = Inf
    end
    # Check if we hit the side
    if norm((point + t_plane * direction)[1:2]) > M.radius
        a = direction[1]^2 + direction[2]^3
        b = 2*(point[1]*direction[1] + point[2]*direction[2])
        c = point[1]^2 + point[2]^2 - M.radius^2
        q = -0.5 * (b - sign(b) * sqrt(b^2 - 4*a*c))
        return sign(q) > 0 ? q/a : c/q
    else # we hit the cap
        return t_plane
    end
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




# RAYRACING
# Generate random free path
random_depth() = -log(rand())

# Generate a random direction weighed with a Henyey-Greenstein
# phase function.
function random_scattering_angle(params::Vector{Float64})
    g = rand() < params[1] ? params[2] : params[3]
    s = 2*rand() - 1
    if g != 0.0
        mu = (1 + g^2 - ((1-g^2)/(1+g*s))^2) / (2g)
    else
        mu = s
    end
    return acos(mu)
end

# Double Henyey-Greenstein function
function double_HG(theta::Real, params::Vector{Float64})
    w = params[1]
    g1 = params[2]
    g2 = params[3]
    ct = cos(theta)
    HG1 = (1 - g1^2) / (1 + g1^2 - 2*g1*ct)^1.5
    HG2 = (1 - g2^2) / (1 + g2^2 - 2*g2*ct)^1.5
    return w*HG1 + (1-w)*HG2
end


# Compute the scattered ray given incoming ray, scattering location, asymmetry parameter
# for the HG phase function, and the single-scattering albedo.
function scattered_ray(ray::Ray, location::Vector{Float64}, phase_params::Vector{Float64}, omega::Float64)
    phi = 2pi * rand()
    theta = random_scattering_angle(phase_params)
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

    t = (ctr*x - str*z)
    new_x = cpr*t - spr*y
    new_y = spr*t + cpr*y
    new_z = str*x + ctr*z
    
    ray.origin = location
    ray.direction[1] = new_x
    ray.direction[2] = new_y
    ray.direction[3] = new_z
    ray.intensity *= omega
    return ray
end

# Function to trace one ray into the medium, starting with intensity 1.0
# incoming from the z-direction. Returns the escaping ray.
function trace_ray(detectors::Array{Detector}, starting_ray::Function, medium::Medium, max_order::Int64, omega::Float64, phase_params::Vector{Float64})
    ray = starting_ray()
    
    # Start the raytracing
    for i = 1:max_order
        tau = random_depth()
        location = ray.origin + tau*ray.direction
        if point_in_medium(medium, location)
            # Peel-off to detector array
            for pixel in detectors
                rho = trace_to_direction(medium, location, pixel.direction)
                scattering_angle = acos(dot(pixel.direction, ray.direction))
                eff = double_HG(scattering_angle, phase_params)
                pixel.intensity[i] += exp(-rho) * ray.intensity * omega * eff
            end
            ray = scattered_ray(ray, location, phase_params, omega)
        else
            return
        end
    end
    return
end


end # Module