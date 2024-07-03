# Geopotential

# Describes a geopotential field and the planetary velocity.
# Instances are callable and return
# the geopotential at a given point in Cartesian coordinates
abstract type Potential end
(potential::Potential)(x) = @inline geopot(potential, x, identity)

# WSW08 geopotential

struct WSW{F} <: Potential
    a::F # Equatorial radius
    g0::F # gravity at r=a when non-rotating (m=eps=0)
    epsilon::F # (a-c)/a
    m::F # non-dimensional squared planetary rotation rate
end
WSW(a, g0, epsilon, m) = WSW(a, g0, epsilon, m)
rotation_rate((; g0, a, m)::WSW) = sqrt(m*g0/a)

function geopot(potential::WSW, xyz, fun=identity)
    (; a, g0, epsilon, m) = potential
    (x,y,z) = fun(xyz)
    R2 = x*x+y*y+z*z
    xi, sinchi2, coschi2 = sqrt(R2)/a, z^2/R2, (x^2+y^2)/R2
    Phi = inv(xi) # order 0
    Phi = Phi - (epsilon-m/2)*Phi^3*(sinchi2-1/3) + (m/2)*coschi2*xi^2 # order 1
    return Phi*(g0*a*a)
end

Earth(T) = WSW(T(1), T(1), T(1/300), T(1/500))

# SGA
struct SGA{F} <: Potential
    a::F # radius
    g0::F # gravity at r=a
    Omega::F # non-dimensional frame rotation rate
end
rotation_rate(potential::SGA) = potential.Omega
SGA(potential) = SGA(potential.a, potential.g0, rotation_rate(potential))

function geopot(potential::SGA, xyz, fun=identity)
    (; a, g0) = potential
    (x,y,z) = fun(xyz)
    return g0*a^2 / sqrt(x*x+y*y+z*z)
end
