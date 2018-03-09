doc"""
apply Hodrick-Prescott filter to `DataVector`.

##### Arguments
- `y::AbstractDataVector` : Data to be detrended
- `λ::Real` : penalty on variation in trend

##### Returns
- `y_df::DataFrame`: `DataFrame` having trend and cyclical data
"""
function hp_filter(y::AbstractDataVector{T}, λ::Real) where T <: Real   
    y_trend, y_cyclical = hp_filter(Vector(y), λ)
    y_df = DataFrame(y_trend = y_trend, y_cyclical = y_cyclical)
    return y_df
end

doc"""
apply Hodrick-Prescott filter to `Vector`.

##### Arguments
- `y::Vector` : `Vector` of values to be detrended
- `λ::Real` : penalty on variation in trend

##### Returns
- `y_trend::Vector`: `Vector` having trend data
- `y_cyclical::Vector`: `Vector` having cyclical data
"""
function hp_filter(y::Vector{T}, λ::Real) where T <: Real
    N = length(y)
    H = hp_filter_matrix(T(λ), N)
    y_trend = H \ y
    y_cyclical = y - y_trend
    return y_trend, y_cyclical
end

doc"""
create a matrix for HP filter

##### Arguments
- `λ::Real` : penalty on variation in trend
- `y::AbstractDataVector` : 

##### Returns
- sparse matrix to be inverted
"""
hp_filter_matrix(λ::Real, N::Integer) =
    spdiagm(-2 => fill(λ, N-2),
            -1 => vcat(-2λ, fill(-4λ, N - 3), -2λ),
             0 => vcat(1 + λ, 1 + 5λ, fill(1 + 6λ, N-4),
                       1 + 5λ, 1 + λ),
             1 => vcat(-2λ, fill(-4λ, N - 3), -2λ),
             2 => fill(λ, N-2))
