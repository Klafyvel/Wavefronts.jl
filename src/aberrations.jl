
"""
Represents an aberration. This uses Zernike polynomials.

See also: [`Piston`](@ref), [`Tilt`](@ref), [`Tip`](@ref), [`Defocus`](@ref), [`ObliqueAstigmatism`](@ref), [`VerticalAstigmatism`](@ref), [`HorizontalComa`](@ref), [`VerticalComa`](@ref), [`ObliqueTrefoil`](@ref), [`VerticalTrefoil`](@ref), [`Spherical`](@ref), [`AddAberration`](@ref), [`Astigmatism`](@ref), [`Coma`](@ref), [`Trefoil `](@ref)
"""
abstract type Aberration end

"""
Piston aberration.
"""
struct Piston <: Aberration
    piston::Real
end
Piston() = Piston(1)

"""
Tilt aberration.
"""
struct Tilt <: Aberration
    tilt::Real
end
Tilt() = Tilt(1)

"""
Tip aberration.
"""
struct Tip <: Aberration
    tip::Real
end
Tip() = Tip(1)

"""
Defocus aberration.
"""
struct Defocus <: Aberration
    defocus::Real
end
Defocus() = Defocus(1)

"""
ObliqueAstigmatism aberration.
"""
struct ObliqueAstigmatism <: Aberration
    oblique::Real
end
ObliqueAstigmatism() = ObliqueAstigmatism(1)

"""
VerticalAstigmatism aberration.
"""
struct VerticalAstigmatism <: Aberration
    vertical::Real
end
VerticalAstigmatism() = VerticalAstigmatism(1)

"""
HorizontalComa aberration.
"""
struct HorizontalComa <: Aberration
    horizontal::Real
end
HorizontalComa() = HorizontalComa(1)

"""
VerticalComa aberration.
"""
struct VerticalComa <: Aberration
    vertical::Real
end
VerticalComa() = VerticalComa(1)

"""
ObliqueTrefoil aberration.
"""
struct ObliqueTrefoil <: Aberration
    oblique::Real
end
ObliqueTrefoil() = ObliqueTrefoil(1)

"""
VerticalTrefoil aberration.
"""
struct VerticalTrefoil <: Aberration
    vertical::Real
end
VerticalTrefoil() = VerticalTrefoil(1)

"""
Spherical aberration.
"""
struct Spherical <: Aberration
    spherical::Real
end
Spherical() = Spherical(1)
abstract type AbstractAddAberration <: Aberration end

"""
The addition of two aberrations. It is created when `+` operator is used on two
[`Aberration`](@ref).
"""
struct AddAberration <: AbstractAddAberration
    left::Aberration
    right::Aberration
end

"""
Astigmatism aberration.
"""
struct Astigmatism <: AbstractAddAberration
    left::ObliqueAstigmatism
    right::VerticalAstigmatism
end
Astigmatism(oblique::Real, vertical::Real) = Astigmatism(ObliqueAstigmatism(oblique), VerticalAstigmatism(vertical))
Astigmatism() = Astigmatism(ObliqueAstigmatism(), VerticalAstigmatism())

"""
Coma aberration.
"""
struct Coma <: AbstractAddAberration
    left::HorizontalComa
    right::VerticalComa
end
Coma(horizontal::Real, vertical::Real) = Coma(HorizontalComa(horizontal), VerticalComa(vertical))
Coma() = Coma(HorizontalComa(), VerticalComa())

"""
Trefoil aberration.
"""
struct Trefoil <: AbstractAddAberration
    left::ObliqueTrefoil
    right::VerticalTrefoil
end
Trefoil(oblique::Real, vertical::Real) = Trefoil(ObliqueTrefoil(oblique), VerticalTrefoil(vertical))
Trefoil() = Trefoil(ObliqueTrefoil(), VerticalTrefoil())

"""
    phase(aberration, coordinates=:cartesian)

Create a function mapping the phase of the `aberration` in the given `coordinates`
system (`:cartesian` or `:polar`).

See also: [`intensity`](@ref), [`amplitude`](@ref)
"""
phase(ab::Aberration, coordinates::Symbol=:cartesian) = phase(convert(Coordinates, coordinates), ab)
phase(t::Coordinates, ab::Aberration) = phase(typeof(t), ab::Aberration)
phase(::Type{Cartesian}, ab::Aberration) = begin
    ϕ = phase(Polar(), ab)
    (x,y) -> begin
        r,θ = cartesian_to_polar(x,y)
        ϕ(r,θ)
    end
end
phase(::Type{Polar}, ab::Piston) = (r,θ) -> ab.piston
phase(::Type{Cartesian}, ab::Piston) = (x,y) -> ab.piston
phase(::Type{Polar}, ab::Tilt) = (r,θ) -> ab.tilt*2r*sin(θ)
phase(::Type{Polar}, ab::Tip) = (r,θ) -> ab.tip*2r*cos(θ)
phase(::Type{Polar}, ab::Defocus) = (r,θ) -> ab.defocus * √3(2r^2-1)
phase(::Type{Polar}, ab::ObliqueAstigmatism) = (r,θ) -> ab.oblique * √6 * r^2 * sin(2θ)
phase(::Type{Polar}, ab::VerticalAstigmatism) = (r,θ) -> ab.vertical * √6 * r^2 * cos(2θ)
phase(::Type{Polar}, ab::HorizontalComa) = (r,θ) -> ab.horizontal * 2*√2*(3r^3-2r)cos(θ)
phase(::Type{Polar}, ab::VerticalComa) = (r,θ) -> ab.vertical * 2*√2*(3r^3-2r)sin(θ)
phase(::Type{Polar}, ab::ObliqueTrefoil) = (r,θ) -> ab.oblique * 2*√2*r^3 * cos(3θ)
phase(::Type{Polar}, ab::VerticalTrefoil) = (r,θ) -> ab.vertical * 2*√2*r^3 * sin(3θ)
phase(::Type{Polar}, ab::Spherical) = (r,θ) -> ab.spherical * √5(6r^4-6r^2+1)
phase(c::Type{Polar}, ab::AbstractAddAberration) = begin
    ϕ1 = phase(c, ab.left)
    ϕ2 = phase(c, ab.right)
    (r,θ) -> ϕ1(r,θ) + ϕ2(r,θ)
end

Base.:+(a::T, b::U) where {T<:Aberration, U<:Aberration} = AddAberration(a,b)

Base.:*(a::T, b::Piston) where {T<:Real} = Piston(a * b.piston)
Base.:*(a::T, b::Tilt) where {T<:Real} = Tilt(a * b.tilt)
Base.:*(a::T, b::Tip) where {T<:Real} = Tip(a * b.tip)
Base.:*(a::T, b::Defocus) where {T<:Real} = Defocus(a * b.defocus)
Base.:*(a::T, b::ObliqueAstigmatism) where {T<:Real} = ObliqueAstigmatism(a * b.oblique)
Base.:*(a::T, b::VerticalAstigmatism) where {T<:Real} = VerticalAstigmatism(a * b.vertical)
Base.:*(a::T, b::HorizontalComa) where {T<:Real} = HorizontalComa(a * b.horizontal)
Base.:*(a::T, b::VerticalComa) where {T<:Real} = VerticalComa(a * b.vertical)
Base.:*(a::T, b::ObliqueTrefoil) where {T<:Real} = ObliqueTrefoil(a * b.oblique)
Base.:*(a::T, b::VerticalTrefoil) where {T<:Real} = VerticalTrefoil(a * b.vertical)
Base.:*(a::T, b::Spherical) where {T<:Real} = Spherical(a * b.spherical)
Base.:*(a::T, b::U) where {T<:Real, U<:AbstractAddAberration} = U(a * b.left, a * b.right)

"""
    intensity(aberration, coordinates=:cartesian; mask=(x,y)->1)

Creates a function mapping the intensity in the given `coordinates` system
(`:cartesian` or `:polar`). This is equivalent of building two `Wavefront`, one
with no aberration and the other with the given `aberration`.

See also: [`phase`](@ref), [`amplitude`](@ref)
"""
intensity(ab::Aberration=Piston(), coordinates::Symbol=:cartesian; mask::Function=((a,b)->1)) = begin
    c = convert(Coordinates, coordinates)
    ϕ = phase(c, ab)
    (a,b) -> abs(1 + mask(a,b)*exp(1im*ϕ(a,b)))^2
end

"""
    project(aberration, data; mask=(x,y)->1)

Project the `data` onto the Zernike polynomial corresponding to the `aberration`
in order to get the fitting [`Aberration`](@ref). A `mask` can be used to fit only
in certain parts of the data (typically usefull if your interference pattern is
only relevant in a circular region). The coordinates are assumed to be normalized
such that (x,y) ∈ [-1;1]².

See also: [`correct`](@ref)
"""
project(ab::Type{<:Aberration}, data::Array{<:Real, 2}; mask::Function=((x,y)->1)) = project(ab(), data, mask=mask)
function project(ab::Aberration, data::Array{<:Real, 2}; mask::Function=((x,y)->1))
    ϕ = phase(Cartesian, ab)
    x,y = range(-1,1,length=size(data)[1]), range(-1,1,length=size(data)[2])
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    mask_values = mask.(X,Y)
    ab_values = ϕ.(X,Y) .* mask_values
    sum(data .* ab_values) / sum(mask_values) * ab
end
function project(ab::AbstractAddAberration, data::Array{<:Real, 2}; mask::Function=((x,y)->1))
    typeof(ab)(project(ab.left, data, mask=mask),project(ab.right, data, mask=mask))
end

"""
    correct(aberration, data; mask=(x,y)->1)

Correct the given aberration in the data by [`project`](@ref)ing the data on the
aberration.

See also: [`project`](@ref)
"""
function correct(ab::Type{<:Aberration}, data::Array{<:Real, 2}; mask::Function=((x,y)->1))
    aberration = project(ab,data,mask=mask)
    ϕ = phase(Cartesian, aberration)
    x,y = range(-1,1,length=size(data)[1]), range(-1,1,length=size(data)[2])
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    mask_values = mask.(X,Y)
    ab_values = ϕ.(X,Y) .* mask_values
    data - ab_values
end

Base.show(a::AddAberration) = begin
    Base.show(a.left)
    Base.show(" + ")
    Base.show(a.right)
end
Base.show(a::Astigmatism) = Base.show("Astigmatism($(a.left.oblique), $(a.right.vertical))")
Base.show(a::Trefoil) = Base.show("Astigmatism($(a.left.oblique), $(a.right.vertical))")
Base.show(a::Coma) = Base.show("Coma($(a.left.horizontal), $(a.right.vertical))")
