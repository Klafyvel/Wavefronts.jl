abstract type Coordinates end
struct Polar <: Coordinates end
struct Cartesian <: Coordinates end

circ(r) = if r < 1 1 else 0 end
circ(::Polar, r,θ) = circ(r)
circ(::Cartesian, x,y) = circ(x^2+y^2)
"""
    circ(x,y)

Circular mask of radius 1.
"""
circ(x,y) = circ(Cartesian(), x, y)

Base.convert(::Type{Coordinates}, s::Symbol) = begin
    if s in [:cartesian :cart :xy]
        Cartesian()
    elseif s in [:polar :pol :rθ]
        Polar()
    else
        error("$s not convertible in coordinates.")
    end
end
