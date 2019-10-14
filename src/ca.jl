# Correspondence Analysis

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


struct CA{T<:Real}

end
