using DSP

function rawphase(
    intensity1::Array{<:Number,2},intensity2::Array{<:Number,2},
    intensity3::Array{<:Number,2},intensity4::Array{<:Number,2}
    )
    atan.(intensity2.-intensity4, intensity1.-intensity3)
end

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
