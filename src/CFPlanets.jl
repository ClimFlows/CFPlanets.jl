module CFPlanets

# API

# the API is WIP. Extending it will be needed to support more approximations,
# e.g. non-traditional / deep atmosphere.

"""
    abstract type ConformalPlanet end

An abstract type for planets such that the two metric factors are equal.
"""
abstract type ConformalPlanet end

"""
    factor = scale_factor(planet::ConformalPlanet, lon, lat) # spherical domain
    factor = scale_factor(planet::ConformalPlanet, x, y)     # logically Cartesian domain

Returns the scale factor at a certain point on `planet``, given by its abstract coordinates.
The scale factor is related to the physical length of a horizontal displacement by:

    δl = scale_factor * sqrt(δlat^2 + cos(lat)^2 δlon^2) # or
    δl = scale_factor * sqrt(δx^2 + δy^2).

For instance if the the spherical-geoid and shallow-atmosphere approximations are made,
the scale factor is the planetary radius, assumed constant.
"""
function scale_factor end

"""
    f = coriolis(planet, lon, lat) # spherical geometry
    f = coriolis(planet, x, y)     # logically Cartesian
Returns the Coriolis parameter corresponding to geometric approximation `planet`.
"""
function coriolis end


#================= Traditional Shallow Spherical geometry ============#

"""
    struct ShallowTradPlanet{F} <: ConformalPlanet
    planet = ShallowTradPlanet(radius, Omega)

Instances of this type describe a spherical-geoid, shallow-atmosphere, traditional
approximation. `radius` is the radius of the planet and `Omega` its rotation rate.
"""
struct ShallowTradPlanet{F} <: ConformalPlanet
    radius::F
    Omega::F
end

@inline scale_factor(planet::ShallowTradPlanet, lon, lat) = planet.radius

@inline coriolis_ST(a, Omega, lat) = 2*Omega*a*a*sin(lat)

@inline coriolis(planet::ShallowTradPlanet, lon, lat) = coriolis_ST.(planet.radius, planet.Omega, lat)

@inline function lonlat_from_cov(ulon, ulat, lon, lat, planet::ShallowTradPlanet)
    invrad = inv(planet.radius)
    return ulon*invrad, ulat*invrad
end

# oblate planets, refactoring needed

include("julia/geopotential.jl")
include("julia/maps.jl")
include("julia/hamiltonian.jl")

#=========================== f-plane in a box ========================#

struct FPlanePlanet{F} <: ConformalPlanet
    dx::F
    f::F
end
@inline scale_factor(planet::FPlanePlanet, x, y) = planet.dx
@inline coriolis(planet::FPlanePlanet, x, y) = planet.f

end # module
