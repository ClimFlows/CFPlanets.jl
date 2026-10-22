"""
    Tank
Abstract type for tanks.
Represents a container with defined dimensions, gravitational acceleration, and the origin `z0` of the geopotential.
"""
abstract type Tank end

"""
    Tank2D
Rectangular 2D tank with length `Lx`, height `Lz`, and gravitational acceleration `g`.
The geopotential is `g (z - z0)`, with `z` the height above the bottom of the tank: `z0` is the height at which it
vanishes (default 0, the bottom).
"""
struct Tank2D{F} <: Tank
    Lx::F
    Lz::F
    g::F
    z0::F
end
Tank2D{F}(Lx, Lz, g) where {F} = Tank2D{F}(Lx, Lz, g, zero(F))
Tank2D(; Lx::F, Lz::F, g::F, z0 = zero(Lz)) where {F} = Tank2D{F}(Lx, Lz, g, z0)

"""
    Tank3D
Rectangular 3D tank with length `Lx`, width `Ly`, height `Lz`, and gravitational acceleration `g`.
The geopotential is `g (z - z0)`, with `z` the height above the bottom of the tank: `z0` is the height at which it
vanishes (default 0, the bottom).
"""
struct Tank3D{F} <: Tank
    Lx::F
    Ly::F
    Lz::F
    g::F
    z0::F
end
Tank3D{F}(Lx, Ly, Lz, g) where {F} = Tank3D{F}(Lx, Ly, Lz, g, zero(F))
Tank3D(; Lx::F, Ly::F, Lz::F, g::F, z0 = zero(Lz)) where {F} = Tank3D{F}(Lx, Ly, Lz, g, z0)
