using DSP

"""
    rawphase(intensity0, intensity1, intensity2, intensity3)

Compute the raw phase from four intensity of interference pattern. The reference
in the pattern of intensity `intensityi` is assumed to be shifted of π/2.

See also: [`unwrapphase`](@ref), [`retrievephase`](@ref)
"""
function rawphase(
    intensity1::Array{<:Number,2},intensity2::Array{<:Number,2},
    intensity3::Array{<:Number,2},intensity4::Array{<:Number,2}
    )
    atan.(intensity2.-intensity4, intensity1.-intensity3)
end

"""
    unwrapphase(wrappedphase)

Unwrap the phase. The method used is inspired by the Standard Experiment Laboratory
script from the MSc in Optics and Photonics at Imperial College London. The input
is typically the output of [`rawphase`](@ref).

See also: [`rawphase`](@ref), [`retrievephase`](@ref)
"""
function unwrapphase(wrappedphase::Array{<:Number,2})
    unwrapped = copy(wrappedphase)
    unwrapped[end÷2+1:end,:] = unwrap(unwrapped[end÷2+1:end, :], dims=1)
    unwrapped = reverse(unwrapped, dims=1)
    unwrapped[end÷2+1:end,:] = unwrap(unwrapped[end÷2+1:end, :], dims=1)
    unwrapped[:, end÷2+1:end] = unwrap(unwrapped[:, end÷2+1:end], dims=2)
    unwrapped = reverse(unwrapped, dims=2)
    unwrapped[:, end÷2+1:end] = unwrap(unwrapped[:, end÷2+1:end], dims=2)
    reverse(reverse(unwrapped, dims=2), dims=1)
end


"""
    retrievephase(intensity0, intensity1, intensity2, intensity3; correct_aberrations, mask)

Basically compose [`rawphase`](@ref) and [`unwrapphase`](@ref). If `correct_aberrations`
is given, tries to [`correct`](@ref) the given aberrations with the given mask.

See also: [`rawphase`](@ref), [`unwrapphase`](@ref), [`correct`](@ref)
"""
function retrievephase(
    intensity1::Array{<:Number,2},intensity2::Array{<:Number,2},
    intensity3::Array{<:Number,2},intensity4::Array{<:Number,2};
    correct_aberrations=[],
    mask=(x,y)->1
    )
    unwrapped = (unwrapphase∘rawphase)(intensity1, intensity2, intensity3, intensity4)
    for aberration in correct_aberrations
        unwrapped = correct(aberration, unwrapped, mask=mask)
    end
    unwrapped
end
