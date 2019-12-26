module Wavefronts


export intensity, field, Wavefront
export Piston, Tilt, Tip, Defocus, Astigmatism, Coma, Trefoil, Spherical, phase, project, correct
export circ
export rawphase, unwrapphase, retrievephase

cartesian_to_polar(x,y) = √(x^2+y^2), atan(x,y)
polar_to_cartesian(r,θ) = r*cos(θ),r*sin(θ)

include("coordinates.jl")
include("aberrations.jl")
include("wavefronts.jl")
include("phasestep.jl")

end # module
