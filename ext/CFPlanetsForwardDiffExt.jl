module CFPlanetsForwardDiffExt

using ForwardDiff: jacobian
import CFPlanets: jac

jac(fwd, (lon,lat,R)) = jacobian(x->[a for a in fwd(x)], [lon, lat, R])

end
