# Correspondence Analysis

# CA type

struct CA{R, S, T<:Real}
    C::CrossTable{R, S, T}  # frequency table
    w_x::Vector{T}          # row weights
    w_y::Vector{T}          # column weights
    M::Matrix{T}            # deviation matrix
    proj::Matrix{T}         # projection matrix
end

function CA(X::Vector{R}, Y::Vector{S})
    C = CrossTable(X, Y)
    n = sum(C)
    C = C / n
    nrows, ncols = size(C)
    w_x = C * ones(eltype(C), ncols)
    w_y = ones(eltype(C), 1, nrows) * C
    M = C - w_x * w_y
    CA(C, w_x, w_y, M)
end

function CA(XY::Matrix{Union{R, S}})
    X = XY[:, 1]
    Y = XY[:, 2]
    CA(X, Y)
end


# Cross tabulation matrix
struct CrossTable{R, S, T<:Integer}
    C::Matrix{T}
    rows::Dict{R, T}
    cols::Dict{S, T}
end

# constructor
function CrossTable(X::Vector{R}, Y::Vector{S}) where {R, S}
    x_n = size(X)
    y_n = size(Y)

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

    CrossTable(table, x_indices, y_indices)
end
