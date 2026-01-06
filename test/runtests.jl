using CFPlanets: WSW, CartesianHamiltonian, SGAHamiltonian,
    OrthogonalHamiltonian, ConformalHamiltonian,
    OblateHamiltonian_I, OblateHamiltonian_II, OblateHamiltonian_III,
    Earth, TD14, TD24, to_curvilinear, to_cartesian, forward,
    metric_factors, orthogonal_factors, conformal_factors,
    Tank2D, Tank3D

using ForwardDiff: jacobian
using StaticArrays
using Test

"""
    xyz, uvw = circular_orbit(planet, R)
Returns an initial position and velocity leading to a circular
orbit of radius R around the planet. Only the spherical part
of geopotential is used for this computation. For non-spherical potential,
the orbit will not be circular.
"""
function circular_orbit(planet, R::F) where {F}
    (; a, g0) = planet
    # acceleration = g0 * (a/R)^2 = omega^2*R
    # velocity = omega*R = a * sqrt(g0/R)
    xyz = normalize(R, (1, 1, 1))
    (ux, uy, uz) = normalize(a * sqrt(g0 / R), vprod(xyz, map(F, (1, -2, 1))))
    return xyz, (ux, uy, uz)
end

function normalize(R, (x, y, z))
    fac = R / sqrt(x^2 + y^2 + z^2)
    return (fac * x, fac * y, fac * z)
end

vprod((x, y, z)::T, (p, q, r)::T) where {T} = T(q * z - r * y, r * x - p * z, p * y - q * x)
vprod((x, y, z)::T, (p, q, r)::T) where {T<:Tuple} =
    q * z - r * y, r * x - p * z, p * y - q * x
dprod((x, y, z), (p, q, r)) = x * p + y * q + z * r

function isconformal(fwd)
    (lon, lat, R) = pi / 5, pi / 4, 1.0
    hlon, hlat, hR = metric_factors(fwd, lon, lat, R)
    @info "Conformality of metric factors" fwd hlon hlat hlon - hlat * cos(lat)^2
end

#===================================== Tests ======================================#

F = Float64
planet = Earth(F)
fwd = TD14(planet)
isconformal(fwd)
fwd = TD24(planet)
isconformal(fwd)

L = 1.
g = 10.

@testset "CFPlanets.jl" begin
    # Write your tests here.

    function check_curvilinear(planet, r)
        x0, p0 = map(SVector{3}, circular_orbit(planet, r))
        xp0 = SVector{6}(x0..., p0...)

        cartesian = CartesianHamiltonian(planet)
        oblateI = OblateHamiltonian_I(planet)
        oblateII = OblateHamiltonian_II(planet)
        oblateIII = OblateHamiltonian_III(planet)
        SGA = SGAHamiltonian(planet)

        h_ref = cartesian(xp0)

        for h in (cartesian, oblateI, oblateII, oblateIII, SGA)
            xi, m = to_curvilinear(h, x0, p0)
            x, p = map(SVector{3}, to_cartesian(h, xi, m))
    #        @info "bwd∘fwd≈identity" x0-x p0-p
            xim = SVector{6}(xi..., m...)
            @info h h(xim) h(xim)-h_ref
            fwd(x) = forward(h, x)
            h1ref, h2ref, h3ref = map(sqrt, metric_factors(fwd, xi...))
            if isa(h, OrthogonalHamiltonian)
                h1, h2, h3, Om, Phi = orthogonal_factors(h, xi[2], xi[3])
                h3 = inv(h3)
                @info "approx scale factors" h1 h2 h3
                @info "error in scale factors" h1-h1ref h2-h2ref h3-h3ref
            elseif isa(h, ConformalHamiltonian)
                h2, h3, Om, Phi = conformal_factors(h, xi[2], xi[3])
                h1, h3 = h2*cos(xi[2]), inv(h3)
                @info "approx scale factors" h1 h2 h3
                @info "error in scale factors" h1-h1ref h2-h2ref h3-h3ref
            end
            println()
        end
    end

    check_curvilinear(WSW(1.0, 1.0, 1e-5, 1e-5), 1.0)

    tank2d = Tank2D(Lx = L, Lz = L, g = g)
    tank3d = Tank3D(Lx = L, Ly = L, Lz = L, g = g)

    @info "Tank2D" tank2d
    @info "Tank3D" tank3d

end
