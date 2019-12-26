abstract type AbstractWavefront end

"""
Represents a wavefront.
"""
struct Wavefront <: AbstractWavefront
    phase::Function
    amplitude::Function
end

"""
Represents the addition of a two wavefronts. It is created when `+` operator is
used on two [`Wavefront`](@ref).
"""
struct AddWavefront <: AbstractWavefront
    left::AbstractWavefront
    right::AbstractWavefront
end
Wavefront(phase::Function, amplitude::Real) = Wavefront(phase, (a,b)->amplitude)
Wavefront(ab::Aberration,coordinates::Symbol=:cartesian) = Wavefront(phase(ab, coordinates), 1)
Wavefront(ab::Aberration, amplitude::Function, coordinates::Symbol=:cartesian) = Wavefront(phase(ab, coordinates), amplitude)
Wavefront() = Wavefront(phase(Piston(0), :cartesian), 1)
Base.:*(num::Number, wavefront::Wavefront) = Wavefront(wavefront.phase, (a,b)->num*wavefront.amplitude(a,b))
Base.:*(wavefront::Wavefront, num::Number) = Wavefront(wavefront.phase, (a,b)->num*wavefront.amplitude(a,b))

Base.:+(wavefront_a::Wavefront, wavefront_b::Wavefront) = AddWavefront(wavefront_a, wavefront_b)

"""
    field(wavefront)

Create a function mapping the field in the coordinates of the wavefront.

See also: [`intensity`](@ref), [`amplitude`](@ref), [`phase`](@ref)
"""
field(wavefront::Wavefront) = (a,b) -> wavefront.amplitude(a,b) * exp(1im*wavefront.phase(a,b))
field(wavefront::AddWavefront) = begin
    field_a = field(wavefront.left)
    field_b = field(wavefront.right)
    (a,b) -> field_a(a,b) + field_b(a,b)
end
field(wavefront::AbstractWavefront, reference::AbstractWavefront) = field(wavefront+reference)

"""
    intensity(wavefront)

Create a function mapping the intensity in the coordinates of the wavefront.

See also: [`field`](@ref), [`amplitude`](@ref), [`phase`](@ref)
"""
intensity(wavefront::AbstractWavefront) = begin
    f = field(wavefront)
    (a,b) -> abs2(f(a,b))
end

"""
    intensity(wavefront, reference)

Create a function mapping the intensity of the interference pattern of the
wavefront with the reference. Both should be in the same coordinate system.

See also: [`field`](@ref), [`amplitude`](@ref), [`phase`](@ref)
"""
intensity(wavefront::AbstractWavefront, reference::AbstractWavefront) = intensity(wavefront+reference)

"""
    phase(wavefront)

Create a function mapping the phase of the wavefront.

See also: [`intensity`](@ref), [`amplitude`](@ref), [`field`](@ref)
"""
phase(wavefront::Wavefront) = wavefront.phase
phase(wavefront::AddWavefront) = begin
    sumf = field(wavefront)
    (a,b) -> angle(sumf(a,b))
end

"""
    amplitude(wavefront)

Create a function mapping the amplitude of the wavefront.

See also: [`intensity`](@ref), [`field`](@ref), [`phase`](@ref)
"""
amplitude(wavefront::Wavefront) = wavefront.amplitude
amplitude(wavefront::AddWavefront) = begin
    sumf = field(wavefront)
    (a,b) -> abs(sumf(a,b))
end
