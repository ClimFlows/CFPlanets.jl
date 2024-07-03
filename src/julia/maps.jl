# Forward and inverse TD14 maps

struct TD14{F}
    planet::WSW{F}
end
struct InvTD14{F}
    planet::WSW{F}
end
Base.inv(fun::TD14) = InvTD14(fun.planet)
Base.inv(fun::InvTD14) = TD14(fun.planet)
(fun::TD14)(xyz) = TD14_fwd(fun.planet, xyz)
(fun::InvTD14)(xyz) = TD14_bwd(fun.planet, xyz)

@fastmath function TD14_bits(planet, lat, xi)
    (; epsilon, m) = planet
    sinlat, coslat = sin(lat), cos(lat)
    A, B = (epsilon-m/2)*inv(xi), m*(xi^2)^2
    dR_E, dR, xi_dchi = A/3 + B/2, A + B/2, (B/3-A)*sinlat*coslat
    return coslat, sinlat, xi, dR_E, dR, xi_dchi
end

@fastmath function TD14_fwd(planet, (lon, lat, R))
    a = planet.a
    coslat, sinlat, xi, dR_E, dR, xi_dchi = TD14_bits(planet, lat, R/a)
    R = xi+dR_E-sinlat^2*dR
    z = R*sinlat + xi_dchi*coslat
    r = R*coslat - xi_dchi*sinlat
    sinlon, coslon = sin(lon), cos(lon)
    return a*r*coslon, a*r*sinlon, a*z
end

@fastmath function TD14_bwd(planet, (x0, y0, z0))
    a = planet.a
    r0, z0 = sqrt(x0^2+y0^2)/a, z0/a
    r1, z1 = r0, z0
    for _ in 1:6
        coslat, sinlat, xi, dR_E, dR, xi_dchi = TD14_bits(planet, atan(z1,r1), sqrt(r1*r1+z1*z1))
        dR = dR_E-sinlat^2*dR
        dr = dR*coslat - xi_dchi*sinlat
        dz = dR*sinlat + xi_dchi*coslat
        r1, z1 = r0-dr, z0-dz # substract O(ε) correction
    end
    return atan(y0,x0), atan(z1,r1), a*sqrt(r1*r1+z1*z1)
end

# Forward and inverse TD24 maps

struct TD24{F}
    planet::WSW{F}
end
struct InvTD24{F}
    planet::WSW{F}
end
Base.inv(fun::TD24) = InvTD24(fun.planet)
Base.inv(fun::InvTD24) = TD24(fun.planet)
(fun::TD24)(xyz) = TD24_fwd(fun.planet, xyz)
(fun::InvTD24)(xyz) = TD24_bwd(fun.planet, xyz)

@fastmath function TD24_fwd(planet, (lon, lat, R))
    (; epsilon, m) = planet
    return TD14_fwd(planet, (lon, lat + (epsilon-5m/6)*cos(lat)*sin(lat), R))
end
@fastmath function TD24_bwd(planet, (x,y,z))
    (; epsilon, m) = planet
    lon, lat14, R = TD14_bwd(planet, (x,y,z))
    # solve lat14 = lat24 + (epsilon-5m/6)*cos(lat24)*sin(lat24)
    lat24 = lat14
    for _ in 1:6
        lat24 = lat14 - (epsilon-5m/6)*cos(lat24)*sin(lat24)
    end
    return lon, lat24, R
end

# metric factors

jac(fwd, (lon,lat,R)) = jacobian(x->[a for a in fwd(x)], [lon, lat, R])

function metric_factors(fwd, lon, lat, R)
    J = jac(fwd, (lon,lat,R))
    g_cov = J'*J # covariant metric tensor
    return g_cov[1,1], g_cov[2,2], g_cov[3,3]
end
