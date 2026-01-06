"""
    Tank
Abstract type for tanks.
Represents a container with defined dimensions and gravitational acceleration.
"""
abstract type Tank end

"""
    Tank2D
Rectangular 2D tank with length `Lx`, height `Lz`, and gravitational acceleration `g`.
"""
struct Tank2D{F} <: Tank
    Lx::F
    Lz::F
    g::F
end
Tank2D(; Lx::F, Lz::F, g::F) where {F} = Tank2D{F}(Lx, Lz, g)

"""
    Tank3D
Rectangular 3D tank with length `Lx`, width `Ly`, height `Lz`, and gravitational acceleration `g`.
"""
struct Tank3D{F} <: Tank
    Lx::F
    Ly::F
    Lz::F
    g::F
end
Tank3D(; Lx::F, Ly::F, Lz::F, g::F) where {F} = Tank3D{F}(Lx, Ly, Lz, g)
