# Correspondence Analysis

# CA type

struct CA{T<:Real}
    mean::Vector{T}         # mean vector
    proj::Matrix{T}         # projection matrix
    prinvars::Vector{T}     # principal variances: of length p
    tprinvar::T             # total principal variance, i.e. sum(prinvars)
    tvar::T                 # total input variance
end

function CA(mean::Vector{T}, proj::Matrix{T}, pvars::Vector{T}, tvar::T) where {T<:Real}
    d, p = size(proj)
    (isempty(mean) || length(mean) == d) ||
        throw(DimensionMismatch("Dimensions of mean and proj are inconsistent."))
    length(pvars) == p ||
        throw(DimensionMismatch("Dimensions of proj and pvars are inconsistent."))
    tpvar = sum(pvars)
    tpvar <= tvar || isapprox(tpvar,tvar) || throw(ArgumentError("principal variance cannot exceed total variance."))
    CA(mean, proj, pvars, tpvar, tvar)
end

function fit(::Type{CA}, X::AbstractMatrix{T};
             method::Symbol=:auto,
             maxoutdim::Int=size(X,1),
             mean=nothing) where {T}

    C, x_indices, y_indices = crosstab(X)
    M, wₘ, wₙ = preprocess_crosstab(C)
    d, n = size(M)

    mv = preprocess_mean(M, mean)

    Z::Matrix{eltype(M)} = centralize(M, mv)

    SVD = svd(Z)
    S = SVD.S
    U = SVD.U

    for i = 1:length(S)
        @inbounds S[i] = abs2(S[i]) / n
    end

    ord = sortperm(S; rev=true)
    var_sum = sum(S)

    CA(mv, U[:, ord], S[ord], var_sum)
end

function transform()
end

mean(M::CA) = M.mean

indim(M::PCA) = size(M.proj, 1)
outdim(M::CA) = size(M.proj, 2)

projection(M::PCA) = M.proj

principalvar(M::PCA, i::Int) = M.prinvars[i]
principalvars(M::PCA) = M.prinvars

transform(M::PCA{T}, x::AbstractVecOrMat{T}) where {T<:Real} = transpose(M.proj) * centralize(x, M.mean)
reconstruct(M::PCA{T}, y::AbstractVecOrMat{T}) where {T<:Real} = decentralize(M.proj * y, M.mean)
