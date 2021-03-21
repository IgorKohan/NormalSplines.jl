function _gram(kernel::RK,
               knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    n_1 = length(knots)
    mat = Matrix{T}(undef, n_1, n_1)
    @inbounds for j = 1:n_1
        for i = 1:j
            mat[i,j] = _rk(kernel, knots[i], knots[j])
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end

function _gram_d1(kernel::RK,
                  knots::Vector{T},
                  d1_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
    n_1 = length(knots)
    n_2 = length(d1_knots)
    n = n_1 + n_2
    mat = Matrix{T}(undef, n, n)
    @inbounds for j = 1:n_1
        for i = 1:j
            mat[i,j] = _rk(kernel, knots[i], knots[j])
            mat[j,i] = mat[i,j]
        end
    end
    n_1_p1 = n_1 + 1
    @inbounds for j = n_1_p1:n
        for i = 1:n_1
            mat[i,j] = _∂rk_∂ξ(kernel, knots[i], d1_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    @inbounds for j = n_1_p1:n
        for i = j:n
            mat[i,j] = _∂²rk_∂η∂ξ(kernel, d1_knots[i-n_1], d1_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end

function _gram_d2(kernel::RK,
                  knots::Vector{T},
                  d2_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
    n_1 = length(knots)
    n_2 = length(d2_knots)
    n = n_1 + n_2
    mat = Matrix{T}(undef, n, n)
    @inbounds for i = 1:n_1
        for j = 1:i
            mat[i,j] = _rk(kernel, knots[i], knots[j])
            mat[j,i] = mat[i,j]
        end
    end
    n_1_p1 = n_1 + 1
    @inbounds for i = 1:n_1
        for j = n_1_p1:n
            mat[i,j] = _∂²rk_∂ξ²(kernel, knots[i], d2_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    @inbounds for i = n_1_p1:n
        for j = i:n
            mat[i,j] = _∂⁴rk_∂η²∂ξ²(kernel, d2_knots[i-n_1], d2_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end

function _gram12(kernel::RK,
               knots::Vector{T},
               d1_knots::Vector{T},
               d2_knots::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
    n_1 = length(knots)
    n_2 = length(d1_knots)
    n_3 = length(d2_knots)
    n = n_1 + n_2 + n_3
    mat = Matrix{T}(undef, n, n)
    @inbounds for j = 1:n_1
        for i = 1:j
            mat[i,j] = _rk(kernel, knots[i], knots[j])
            mat[j,i] = mat[i,j]
        end
    end
    n_1_p1 = n_1 + 1
    ne = n_1 + n_2
    @inbounds for j = n_1_p1:ne
        for i = 1:n_1
            mat[i,j] = _∂rk_∂ξ(kernel, knots[i], d1_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    ns = n_1 + n_2 + 1
    @inbounds for j = ns:n
        for i = 1:n_1
            mat[i,j] = _∂²rk_∂ξ²(kernel, knots[i], d2_knots[j-ns+1])
            mat[j,i] = mat[i,j]
        end
    end
    ne = n_1 + n_2
    @inbounds for j = n_1_p1:ne
        for i = j:ne
            mat[i,j] = _∂²rk_∂η∂ξ(kernel, d1_knots[i-n_1], d1_knots[j-n_1])
            mat[j,i] = mat[i,j]
        end
    end
    ns = n_1 + n_2 + 1
    ne = n_1 + n_2
    @inbounds for j = ns:n
        for i = n_1_p1:ne
            mat[i,j] = _∂³rk_∂η∂ξ²(kernel, d1_knots[i-n_1], d2_knots[j-ns+1])
            mat[j,i] = mat[i,j]
        end
    end
    ns = n_1 + n_2 + 1
    @inbounds for j = ns:n
        for i = j:n
            mat[i,j] = _∂⁴rk_∂η²∂ξ²(kernel, d2_knots[i-ns+1], d2_knots[j-ns+1])
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end
