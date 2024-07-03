#================== Hamiltonians ====================#

# `h::Hamiltonian` is a callable that yields H((x,y,z,p,q,r))
# its gradient will be taken when solving the ODE of motion.
abstract type Hamiltonian end

# a Hamiltonian also knows how to map to/from Cartesian coordinates
function forward end # maps from curvilinear to Cartesian
function backward end # maps from Cartesian to curvilinear

#===== Hamiltonian in Cartesian coordinates, in the rotating frame =====#

# "planet" knows about geopotential Φ(x,y,z) and rotation rate Ω
struct CartesianHamiltonian{P} <: Hamiltonian
    planet::P
end
# returns Hamiltonian (scalar) given position-momentum (x,y,z,p,q,r)
# in Cartesian coordinates (p,q,r) = absolute velocity
# => substract planetary velocity to get relative velocity => kinetic energy
function (h::CartesianHamiltonian)((x, y, z, p, q, r))
    Omega = rotation_rate(h.planet)
    u, v, w = p + Omega * y, q - Omega * x, r # substract planetary velocity
    return (u^2 + v^2 + w^2) / 2 - h.planet((x, y, z)) # NB : geopotential has unusual sign convention
end

forward(::CartesianHamiltonian, x) = x
backward(::CartesianHamiltonian, x) = x

#============= Hamiltonian in orthogonal geopotential coordinates ===========#

abstract type OrthogonalHamiltonian <: Hamiltonian end
"""
    h₁, h₂, h³, R¹, Φ = orthognoal_factors(h::OrthogonalHamiltonian, ξ¹, ξ², ξ³)
"""
function orthogonal_factors end

function (h::OrthogonalHamiltonian)((λ, ϕ, ξ³, m₁, m₂, m₃))
    h₁, h₂, h³, R¹, Φ = orthogonal_factors(h, ϕ, ξ³)
    u = m₁ / h₁ - R¹*h₁
    v = m₂ / h₂
    w = m₃ * h³
    return (u^2+v^2+w^2)/2 - Φ # kinetic + potential
end

#============= oblate approximation I: O(ε²) errors ===========#

struct OblateHamiltonian_I{F} <: OrthogonalHamiltonian
    planet::WSW{F}
    Omega::F
end
OblateHamiltonian_I(p::WSW) = OblateHamiltonian_I(p, rotation_rate(p))
forward(h::OblateHamiltonian_I, lonlatxi) = TD14_fwd(h.planet, lonlatxi)
backward(h::OblateHamiltonian_I, xyz) = TD14_bwd(h.planet, xyz)

function orthogonal_factors(h::OblateHamiltonian_I, phi, xi)
    (; epsilon, m, g0, a) = h.planet
    sinphi, cosphi = sincos(phi)
    sinphi2, cos2phi = sinphi^2, cosphi^2 - sinphi^2
    # ξ is a pseudo-radius (dimension = length, like a)
    A, B = (epsilon - m / 2) * (a / xi)^2, m * (xi / a)^3
    # horizontal metric factors, dimension = length
    R_E, Delta_R, X = xi * (1 + A / 3 + B / 2), xi * (A + B / 2), xi * (A - B / 3)
    h_lambda = cosphi * (R_E - Delta_R * sinphi2 + X * sinphi2)
    h_phi = R_E - Delta_R * sinphi2 - X * cos2phi
    # h_xi = ξ²g = contravariant vertical metric factor, close to 1
    h_xi = 1 + ((A / 3 - 2B) + (2B - A) * sinphi2)
    return h_lambda, h_phi, h_xi, h.Omega, inv(xi) * (g0 * a)
end

#============= oblate approximation II: TD24 mapping ============#

struct OblateHamiltonian_II{F} <: OrthogonalHamiltonian
    planet::WSW{F}
    Omega::F
end
OblateHamiltonian_II(p::WSW) = OblateHamiltonian_II(p, rotation_rate(p))
forward(h::OblateHamiltonian_II, lonlatxi) = TD24_fwd(h.planet, lonlatxi)
backward(h::OblateHamiltonian_II, xyz) = TD24_bwd(h.planet, xyz)

function orthogonal_factors(h::OblateHamiltonian_II, phi, xi)
    (; epsilon, m, g0, a) = h.planet
    sinphi, cosphi = sincos(phi)
    sinphi2, cos2phi = sinphi^2, cosphi^2 - sinphi^2
    # ξ is a pseudo-radius (dimension = length, like a)
    A, B = (epsilon - m / 2) * (a / xi)^2, m * (xi / a)^3
    # horizontal metric factors, dimension = length
    R_E, Delta_R = xi * (1 + A / 3 + B / 2), xi * (A + B / 2)
    Y = xi * (A - B / 3 - (epsilon-(5m/6)))
    h_lambda = cosphi * (R_E - Delta_R * sinphi2 + Y * sinphi2)
    h_phi = R_E - Delta_R * sinphi2 - Y * cos2phi
    # h_xi = ξ²g = contravariant vertical metric factor, close to 1
    h_xi = 1 + ((A / 3 - 2B) + (2B - A) * sinphi2)
    return h_lambda, h_phi, h_xi, h.Omega, inv(xi) * (g0 * a)
end

#============= Hamiltonian in conformal geopotential coordinates ===========#

abstract type ConformalHamiltonian <: Hamiltonian end
"""
    h₂, h³, R¹, Φ = conformal_factors(h::OrthogonalHamiltonian, ξ¹, ξ², ξ³)
"""
function conformal_factors end

function (h::ConformalHamiltonian)((λ, ϕ, ξ³, m₁, m₂, m₃))
    h₂, h³, R¹, Φ = conformal_factors(h, ϕ, ξ³)
    h₁ = h₂ * cos(ϕ)
    u = m₁ / h₁ - R¹*h₁
    v = m₂ / h₂
    w = m₃ * h³
    return (u^2+v^2+w^2)/2 - Φ # kinetic + potential
end

#============= oblate approximation III: simplified oblate deep-atmosphere ============#

struct OblateHamiltonian_III{F} <: ConformalHamiltonian
    planet::WSW{F}
    Omega::F
    xi_ref::F
end
OblateHamiltonian_III(p::WSW) = OblateHamiltonian_III(p, rotation_rate(p), xi_ref(p))
xi_ref(p::WSW) = 1 - (p.epsilon + p.m) / 3

forward(h::OblateHamiltonian_III, (lon, lat, zeta)) =
    TD24_fwd(h.planet, (lon, lat, zeta * h.xi_ref))

function backward(h::OblateHamiltonian_III, xyz)
    (lon, lat, xi) = TD24_bwd(h.planet, xyz)
    return (lon, lat, xi / h.xi_ref)
end

function conformal_factors(h::OblateHamiltonian_III, phi, zeta)
    (; planet, Omega, xi_ref) = h
    (; epsilon, m) = planet
    sinphi2 = sin(phi)^2
    h = zeta * (1 - epsilon * sinphi2)
    g_ref = (1 + (epsilon - 3m / 2)) * (1 + (5m / 2 - epsilon) * sinphi2)
    return h, xi_ref*g_ref, Omega, inv(zeta * xi_ref)
end

#============= SGA ===========#

struct SGAHamiltonian{F} <: ConformalHamiltonian
    a::F
    g0::F
    Omega::F
    planet::WSW{F}
end
SGAHamiltonian(p::WSW) = SGAHamiltonian(p.a, p.g0, rotation_rate(p), p)

# we could also use TD24
forward(h::SGAHamiltonian, lonlatxi) = TD14_fwd(h.planet, lonlatxi)
backward(h::SGAHamiltonian, xyz) = TD14_bwd(h.planet, xyz)

conformal_factors(h::SGAHamiltonian, _, R) = R, one(R), h.Omega, (h.g0 * h.a^2) / R

#================== Hamiltonian dynamics ====================#

struct HamiltonDynamics{H<:Hamiltonian}
    hamiltonian::H
end

#=
function (dyn::HamiltonDynamics)(dxp::MVector{6}, (x, y, z, p, q, r)::MVector{6}, _, t)
    xp = SVector(x, y, z, p, q, r)
    Hx, Hy, Hz, Hp, Hq, Hr = gradient(dyn.hamiltonian, xp)
    dxp[1], dxp[2], dxp[3] = Hp, Hq, Hr
    dxp[4], dxp[5], dxp[6] = -Hx, -Hy, -Hz
end
=#
