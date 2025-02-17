
# print arrays in pretty way

function printarr(io::IO, a::AbstractArray)
    Base.with_output_limit(()->Base.showarray(io, a, header=false, repr=false))
end

printvec(io::IO, a::AbstractVector) = printarr(io, a')

printarrln(io::IO, a::AbstractArray) = (printarr(io, a); println(io))
printvecln(io::IO, a::AbstractVector) = (printvec(io, a); println(io))

# centralize

centralize(x::AbstractVector, m::AbstractVector) = (isempty(m) ? x : x - m)
centralize(x::AbstractMatrix, m::AbstractVector) = (isempty(m) ? x : x .- m)

decentralize(x::AbstractVector, m::AbstractVector) = (isempty(m) ? x : x + m)
decentralize(x::AbstractMatrix, m::AbstractVector) = (isempty(m) ? x : x .+ m)

# get a full mean vector

fullmean(d::Int, mv::Vector{T}) where T = (isempty(mv) ? zeros(T, d) : mv)::Vector{T}

preprocess_mean(X::AbstractMatrix{T}, m) where T<:Real =
    (m == nothing ? vec(mean(X, dims=2)) : m == 0 ? T[] :  m)::Vector{T}

# choose the first k values and columns
#
# S must have fields: values & vectors

function extract_kv(fac::Factorization{T}, ord::AbstractVector{Int}, k::Int) where T
    si = ord[1:k]
    vals = fac.values[si]::Vector{T}
    vecs = fac.vectors[:, si]::Matrix{T}
    return (vals, vecs)
end


# symmmetrize a matrix

function symmetrize!(A::Matrix)
    n = size(A, 1)
    @assert size(A, 2) == n
    for j = 1:n
        for i = 1:j-1
            @inbounds A[i,j] = A[j,i]
        end
        for i = j+1:n
            @inbounds A[i,j] = middle(A[i,j], A[j,i])
        end
    end
    return A
end

# percolumn dot

function coldot(X::AbstractMatrix{T}, Y::AbstractMatrix{T}) where T<:Real
    m = size(X, 1)
    n = size(X, 2)
    @assert size(Y) == (m, n)
    R = zeros(T, n)
    for j = 1:n
        R[j] = dot(view(X,:,j), view(Y,:,j))
    end
    return R
end

# qnormalize!

function qnormalize!(X, C)
    # normalize each column of X (say x), such that x'Cx = 1
    m = size(X, 1)
    n = size(X, 2)
    CX = C * X
    for j = 1:n
        x = view(X,:,j)
        cx = view(CX,:,j)
        rmul!(x, inv(sqrt(dot(x, cx))))
    end
    return X
end

# add_diag!

function add_diag!(A::AbstractMatrix, v::Real)
    # add v to diagonal of A
    m = size(A, 1)
    n = size(A, 2)
    @assert m == n
    if v != zero(v)
        for i = 1:n
            @inbounds A[i,i] += v
        end
    end
    return A
end

# regularize a symmetric matrix
function regularize_symmat!(A::Matrix{T}, lambda::Real) where T<:Real
    if lambda > 0
        emax = eigmax(Symmetric(A))
        add_diag!(A, emax * lambda)
    end
    return A
end

"""
    calcscattermat([covestimator::CovarianceEstimator], Z::DenseMatrix)

Calculate the scatter matrix of centered data `Z` based on a covariance
matrix calculated using covariance estimator `covestimator` (by default,
`SimpleCovariance()`).
"""
function calcscattermat(covestimator::CovarianceEstimator, Z::DenseMatrix{T}) where T<:Real
    return cov(covestimator, Z; dims=2, mean=zeros(T, size(Z, 1)))*size(Z, 2)
end

function calcscattermat(Z::DenseMatrix)
    return calcscattermat(SimpleCovariance(), Z)
end

# calculate pairwise kernel
function pairwise!(K::AbstractVecOrMat{T}, kernel::Function,
                   X::AbstractVecOrMat{T}, Y::AbstractVecOrMat{T}) where T<:AbstractFloat
    n = size(X, 2)
    m = size(Y, 2)
    for j = 1:m
        aj = view(Y, :, j)
        for i in j:n
            @inbounds K[i, j] = kernel(view(X, :, i), aj)[]
        end
        j <= n && for i in 1:(j - 1)
            @inbounds K[i, j] = K[j, i]   # leveraging the symmetry
        end
    end
    K
end

pairwise!(K::AbstractVecOrMat{T}, kernel::Function, X::AbstractVecOrMat{T}) where T<:AbstractFloat =
    pairwise!(K, kernel, X, X)

function pairwise(kernel::Function, X::AbstractVecOrMat{T}, Y::AbstractVecOrMat{T}) where T<:AbstractFloat
    n = size(X, 2)
    m = size(Y, 2)
    K = similar(X, n, m)
    pairwise!(K, kernel, X, Y)
end

pairwise(kernel::Function, X::AbstractVecOrMat{T}) where T<:AbstractFloat =
    pairwise(kernel, X, X)

# calculate cross tabulation of 2 variables
function crosstab(X::Matrix{R}) where R
    nrows, ncols = size(X)
    @assert ncols == 2 "Use a burt table for more than 2 variables"

    crosstab(X[:, 1], X[:, 2])
end

function crosstab(X::Vector{R}, Y::Vector{S}) where {R, S}
    x_n = length(X)
    y_n = length(Y)

    @assert x_n == y_n "Both vectors must contain the same number of elements"

    x_values = unique(X)
    y_values = unique(Y)

    x_values_n = length(x_values)
    y_values_n = length(y_values)

    x_indices = Dict{R, Int}(zip(x_values, 1:x_values_n))
    y_indices = Dict{S, Int}(zip(y_values, 1:y_values_n))

    table = zeros(Int, x_values_n, y_values_n)

    for pair in zip(X, Y)
        x_index = x_indices[pair[1]]
        y_index = y_indices[pair[2]]
        table[x_index, y_index] += 1
    end

    return table, x_indices, y_indices
end

# Calculate χ² statistic
function preprocess_crosstab(C::Matrix{T}) where {T<:Int}
    n = sum(C)
    # Frequency matrix
    F = C / n
    nrows, ncols = size(C)

    # joint probability distribution
    wₘ = F * ones(eltype(F), ncols)
    wₙ = ones(eltype(F), 1, nrows) * F

    # χ² statistic
    M = F - wₘ * wₙ

    return M, wₘ, wₙ
end
