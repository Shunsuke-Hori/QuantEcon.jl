hp_filter(y::AbstractVector{T}, λ::Real) where T <: Real = hp_filter(Vector(y), λ)

function hp_filter(y::Vector{T}, λ::Real) where T <: Real
    N = length(y)
    H = hp_filter_matrix(T(λ), N)
    y_trend = H \ y
    y_cyclical = y - y_trend
    return y_trend, y_cyclical
end

hp_filter_matrix(λ::Real, N::Integer) =
    spdiagm(-2 => fill(λ, N-2),
            -1 => vcat(-2λ, fill(-4λ, N - 3), -2λ),
             0 => vcat(1 + λ, 1 + 5λ, fill(1 + 6λ, N-4),
                       1 + 5λ, 1 + λ),
             1 => vcat(-2λ, fill(-4λ, N - 3), -2λ),
             2 => fill(λ, N-2))