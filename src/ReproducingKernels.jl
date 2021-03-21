
@doc raw"
`struct RK_W1 <: ReproducingKernel_1`

Define a type of reproducing kernel of Sobolev space ``W^2_1 [0,1]``:
```math
V(\eta, \xi) =
  \begin{cases}
       1 + \eta \,  , &  0 \le \eta \le \xi \le 1 \\
       1 + \xi \,  , &  0 \le \xi \le \eta \le 1
    \end{cases}
```
"
struct RK_W1 <: ReproducingKernel_1
   RK_W1() = new()
   RK_W1(p::Int64) = new()
   RK_W1(p::AbstractFloat) = new()
end

@doc raw"
`struct RK_W2 <: ReproducingKernel_2`

Define a type of reproducing kernel of Sobolev space ``W^2_2 [0,1]``:
```math
V(\eta, \xi) =
  \begin{cases}
       1 + (\eta + \eta^2 / 2)\xi - \eta^3 / 6  \,  , &  0 \le \eta \le \xi \le 1  \\
       1 + (\xi + \xi^2 / 2)\eta - \xi^3 / 6 \,  , &  0 \le \xi \le \eta \le 1
    \end{cases}
```
"
struct RK_W2 <: ReproducingKernel_2
   RK_W2() = new()
   RK_W2(p::Int64) = new()
   RK_W2(p::AbstractFloat) = new()
end

@doc raw"
`struct RK_W3 <: ReproducingKernel_3`

Define a type of reproducing kernel of Sobolev space ``W^2_3 [0,1]``:
```math
V(\eta, \xi) =
  \begin{cases}
       \sum_{i=0}^2 \frac{\xi^i}{!i} \left ( \frac{\eta^i}{!i} + (-1)^{i} \frac {\eta^{5 - i}}{(5 - i)!}  \right )  \,  , &  0 \le \eta \le \xi \le 1
        \\
       \sum_{i=0}^2 \frac{\eta^i}{!i} \left ( \frac{\xi^i}{!i} + (-1)^{i} \frac {\xi^{5 - i}}{(5 - i)!}  \right ) \,  , &  0 \le \xi \le \eta \le 1   \end{cases}
```
"
struct RK_W3 <: ReproducingKernel_3
   RK_W3() = new()
   RK_W3(p::Int64) = new()
   RK_W3(p::AbstractFloat) = new()
end

#
@doc raw"
`struct RK_H1{T <: AbstractFloat} <: ReproducingKernel_1`

Define a type of reproducing kernel of Bessel Potential space ``H^1_ε (R)``:
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|) \, .
```
# Arguments
- no argument: means that parameter `ε` will be estimated in course of
               interpolation procedure execution.
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it must be greater than zero.
"
struct RK_H1{T <: AbstractFloat} <: ReproducingKernel_1
     ε::T
     RK_H1() = new{Float64}(0.0)
     function RK_H1(ε::T) where T <: AbstractFloat
        if ε <= 0
          throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H1(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@doc raw"
`struct RK_H2{T <: AbstractFloat} <: ReproducingKernel_2`

Define a type of reproducing kernel of Bessel Potential space ``H^2_ε (R)``:
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|)
             (1 + \varepsilon |\xi  - \eta|) \, .
```
# Arguments
- no argument: means that parameter `ε` will be estimated in course of
               interpolation procedure execution.
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it must be greater than zero.
"
struct RK_H2{T <: AbstractFloat} <: ReproducingKernel_2
     ε::T
     RK_H2() = new{Float64}(0.0)
     function RK_H2(ε::T) where T <: AbstractFloat
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H2(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@doc raw"
`struct RK_H3{T <: AbstractFloat} <: ReproducingKernel_3`

Define a type of reproducing kernel of Bessel Potential space ``H^3_ε (R)``:
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|)
             (3 + 3\varepsilon |\xi  - \eta| + \varepsilon ^2 |\xi - \eta| ^2 ) \, .
```
# Arguments
- no argument: means that parameter `ε` will be estimated in course of
               interpolation procedure execution.
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it must be greater than zero.
"
struct RK_H3{T <: AbstractFloat} <: ReproducingKernel_3
     ε::T
     RK_H3() = new{Float64}(0.0)
     function RK_H3(ε::T) where T <: AbstractFloat
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H3(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@inline function _rk(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
   value::T = 0.0
   defined::Bool  = false
   if η <= ξ
      if isa(kernel, RK_H3)
         defined = true
         x = kernel.ε * (ξ - η)
         value = (3 + x * (3 + x)) * exp(-x)
      end
      if isa(kernel, RK_H2)
         defined = true
         x = kernel.ε * (ξ - η)
         value = (1 + x) * exp(-x)
      end
      if isa(kernel, RK_H1)
         defined = true
         x = kernel.ε * (ξ - η)
         value = exp(-x)
      end
      if isa(kernel, RK_W3)
         defined = true
         value = 1 + η * (ξ + η * (η^2 * (η - 5 * ξ) + 10 * ξ^2 * (η + 3))/120)
      end
      if isa(kernel, RK_W2)
         defined = true
         value = 1 + (1 + η/2) * η * ξ - η^3/6
      end
      if isa(kernel, RK_W1)
         defined = true
         value = 1 + η
      end
   else
      value = _rk(kernel, ξ, η)
      return value
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂rk_∂ξ(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = kernel.ε * (ξ - η)
      value = -kernel.ε * x * (1 + abs(x)) * exp(-abs(x))
   end
   if isa(kernel, RK_H2)
      defined = true
      x = kernel.ε * (ξ - η)
      value = -kernel.ε * x * exp(-abs(x))
   end
   if isa(kernel, RK_W3)
       defined = true
       if η <= ξ
          value = η + η * (4 * η * ξ * (η + 3) - η^3)/24
       else
          value = (ξ^4 - 4 * η * (ξ^3 - 6) + 6 * ξ * η^2 * (ξ + 2))/24
       end
   end
   if isa(kernel, RK_W2)
      defined = true
      if η <= ξ
         value = (1 + η/2) * η
      else
         value = ξ * (η - ξ/2) + η
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂²rk_∂η∂ξ(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
   value::T = 0.0
   defined::Bool  = false
   if η <= ξ
      if isa(kernel, RK_H3)
         defined = true
         x = abs(kernel.ε * (ξ - η))
         value = -(kernel.ε^2) * (x * (x - 1) - 1) * exp(-x)
      end
      if isa(kernel, RK_H2)
         defined = true
         x = abs(kernel.ε * (ξ - η))
         value = kernel.ε^2 * (1 - x) * exp(-x)
      end
      if isa(kernel, RK_W3)
         defined = true
         value = -(η^3)/6 + η * ξ * (η + 2)/2 + 1
      end
      if isa(kernel, RK_W2)
         defined = true
         value = 1 + η
      end
   else
      value = _∂²rk_∂η∂ξ(kernel, ξ, η)
      return value
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂²rk_∂ξ²(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = abs(kernel.ε * (ξ - η))
      value = kernel.ε^2 * (x * (x - 1) - 1) * exp(-x)
   end
   if isa(kernel, RK_W3)
      defined = true
      if η <= ξ
         value = η^2 * (η + 3)/6
      else
         value = (ξ^2 * (ξ - 3 * η) + 3 * η^2 * (ξ + 1))/6
      end
   end

   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂³rk_∂η∂ξ²(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = kernel.ε * (ξ - η)
      value = kernel.ε^3 * (x * (abs(x) - 3)) * exp(-abs(x))
   end
   if isa(kernel, RK_W3)
      defined = true
      if η <= ξ
          value = η * (1 + η/2)
      else
         value = ξ * (η - ξ/2) + η
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂⁴rk_∂η²∂ξ²(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = abs(kernel.ε * (ξ - η))
      value = kernel.ε^4 * (x * (x - 5) + 3) * exp(-x)
   end
   if isa(kernel, RK_W3)
      defined = true
      value = η + 1
      if η <= ξ
         value = η + 1
      else
         value = ξ + 1
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂rk_∂η(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_2}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = kernel.ε * (ξ - η)
      value = kernel.ε * x * (1 + abs(x)) * exp(-abs(x))
   end
   if isa(kernel, RK_H2)
      defined = true
      x = kernel.ε * (ξ - η)
      value = kernel.ε * x * exp(-abs(x))
   end
   if isa(kernel, RK_W3)
       defined = true
       if η <= ξ
          value = (η * (η^3 + 6 * ξ^2 * (η + 2)) - 4 * ξ * (η^3 - 6))/24
       else
          value = ξ + ξ * (4 * η * ξ * (ξ + 3) - ξ^3)/24
       end
   end
   if isa(kernel, RK_W2)
      defined = true
      if η <= ξ
         value = η * (ξ - η/2) + ξ
      else
         value = (1 + ξ/2) * ξ
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂²rk_∂η²(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = abs(kernel.ε * (ξ - η))
      value = kernel.ε^2 * (x * (x - 1) - 1) * exp(-x)
   end
   if isa(kernel, RK_W3)
      defined = true
      if η <= ξ
         value = (η^2 * (η - 3 * ξ) + 3 * ξ^2 * (η + 1))/6
      else
         value = ξ^2 * (ξ + 3)/6
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end

@inline function _∂³rk_∂η²∂ξ(kernel::RK, η::T, ξ::T) where {T <: AbstractFloat, RK <: ReproducingKernel_3}
   value::T = 0.0
   defined::Bool  = false
   if isa(kernel, RK_H3)
      defined = true
      x = kernel.ε * (ξ - η)
      value = -kernel.ε^3 * (x * (abs(x) - 3)) * exp(-abs(x))
   end
   if isa(kernel, RK_W3)
      defined = true
      if η <= ξ
          value = η * (ξ - η/2) + ξ
      else
          value = ξ * (1 + ξ/2)
      end
   end
   if !defined
      Throw(ArgumentError("kernel: Incorrect parameter type."))
   end
   return value
end
