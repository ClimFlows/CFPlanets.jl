# CFPlanets

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ClimFlows.github.io/CFPlanets.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ClimFlows.github.io/CFPlanets.jl/dev/)
[![Build Status](https://github.com/ClimFlows/CFPlanets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClimFlows/CFPlanets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ClimFlows/CFPlanets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ClimFlows/CFPlanets.jl)

The aim of the `CFPlanets` is to provide a single source of truth for common and less common geometric-dynamical approximations, i.e:
* common: traditional shallow-atmosphere approximation, f-plane
* less common: non-traditional f-plane and shallow-atmosphere, deep-atmosphere, ...

## Installation

`CFPlanets` is registered in the ClimFlows [repository](https://github.com/ClimFlows/JuliaRegistry). Follow instructions there, then:
```julia
] add CFPlanets
```

## Status

The general idea behind `CFPlanets` is to formulate geophysical models in an abstract coordinate system devoid of metric and dynamical information/hypotheses (e.g. the unit sphere parameterized by latitude and longitude). Physical equations written in such coordinates then explicitly incorporate metric and Coriolis information reflecting the specific geometric-dynamical approximation being made. `CFPlanets` is here to compute the required information (e.g. metric factors).

The API is not stabilized yet and subject to change.
