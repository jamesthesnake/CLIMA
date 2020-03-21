using MPI
using Test
using LinearAlgebra
using Random
using GPUifyLoops, StaticArrays
using CLIMA
using CLIMA.LinearSolvers
using CLIMA.ColumnwisePreconditionedConjugateGradientSolver
using CuArrays

using Revise

CLIMA.init()
const ArrayType = CLIMA.array_type()
#const device = ArrayType == Array ? CPU() : CUDA()
const device = CUDA()

n = 100
T  = Float64
A = rand( n, n)
const scale = 2.0
const ϵ = 0.1
# Matrix 1
A = A' * A .* ϵ + scale*I

# Matrix 2
A = CuArray(Diagonal(collect(1:n) * 1.0))
#=
positive_definite = minimum(eigvals(A)) > eps(1.0)
if positive_definite
  println("the matrix is positive definite")
else
  println("The matrix is not positive definite")
end
=#

b = CuArray(ones(n)*1.0)
mulbyA!(y, x) = (y .= A * x)

tol = sqrt(eps(T))
method(b,tol) = ColumnwisePreconditionedConjugateGradient(b, max_iter=n)
linearsolver = method(b, tol)

x = CuArray(ones(n)*1.0)
x0 = copy(x)
println("------------------------")
println("------------------------")
iters = linearsolve!(mulbyA!, linearsolver, x, b; max_iters=Inf)
exact = A \ b
x0 = copy(x)
println("------------------------")
println("------------------------")
println("the norm of the error is ")
println(norm(x - exact)/norm(exact))
println("the relative norm of the residual is ")
println(norm(A*x - b)/norm(b))
println("The number of iterations is ")
println(iters)

###

# A few structures for defining more complex linear operators
function entry(i1, i2, i3, i4, i5, i6; FT = Float64)
  val = FT((i1*i2)*(1 + i3 + i5))
  return val
end

# A few structures for defining more complex linear operators
#=
function entry(i1, i2, i3, i4, i5, i6; FT = Float64)
  val = FT(1.0)
  return val
end
=#
# A few structures for defining more complex linear operators
function inverse_entry(i1, i2, i3, i4, i5, i6; FT = Float64)
  val = 1.0 / entry(i1, i2, i3, i4, i5, i6; FT = FT)
  return val
end

"""
more_complicated_linear_operator!(y,x)

# Description

- Testing columnwise solve and reduction for a simple operator

- Acts on a vector "a" with size(a) = (n_i, n_j, n_k, n_s, n_ev, n_eh, n_r)

- This is a diagonal linear operator

# Arguments
- `y`: (array). OVERWRITTEN. output of linear operator acting on x
- `x`: (array). thing that the linear operator acts on
"""
function more_complicated_linear_operator!(y,x)
  tup = size(y)
  for i6 in 1:tup[6]
    for i5 in 1:tup[5]
      for i4 in 1:tup[4]
        for i3 in 1:tup[3]
          for i2 in 1:tup[2]
            for i1 in 1:tup[1]
              y[i1, i2, i3, i4, i5, i6] = entry(i1, i2, i3, i4, i5, i6) * x[i1, i2, i3, i4, i5, i6]
            end
          end
        end
      end
    end
  end
end

"""
more_complicated_inverse_linear_operator!(y,x)

# Description

- Testing columnwise solve and reduction for a simple operator

- Acts on a vector "a" with size(a) = (n_i, n_j, n_k, n_s, n_ev, n_eh, n_r)

- This is a diagonal linear operator

# Arguments
- `y`: (array). OVERWRITTEN. output of linear operator acting on x
- `x`: (array). thing that the linear operator acts on
"""
function more_complicated_inverse_linear_operator!(y,x)
  tup = size(y)
  for i6 in 1:tup[6]
    for i5 in 1:tup[5]
      for i4 in 1:tup[4]
        for i3 in 1:tup[3]
          for i2 in 1:tup[2]
            for i1 in 1:tup[1]
              y[i1, i2, i3, i4, i5, i6] = inverse_entry(i1, i2, i3, i4, i5, i6) * x[i1, i2, i3, i4, i5, i6]
            end
          end
        end
      end
    end
  end
end
tup = (3,4,5, 1, 2, 1)
y = randn(tup)
x = copy(y)
z = copy(x)

more_complicated_linear_operator!(y,x)
more_complicated_inverse_linear_operator!(x, y)
norm(z .- x)


method(x,tol) = ColumnwisePreconditionedConjugateGradient(x, max_iter= tup[3]*tup[5], dims = (3,5))
linearsolver = method(x, tol)

y = copy(x)
iters = linearsolve!(more_complicated_linear_operator!, linearsolver, x, y; max_iters=Inf)
exact = copy(x)
more_complicated_inverse_linear_operator!(exact, y)
println("the norm of the error is ")
println(norm(x - exact)/norm(exact))
println("The number of iterations is ")
println(iters)


###
# Test with something more akin to CliMA structure

"""
closure_even_more_complicated_linear_operator!(A)

# Description
- Closure for even_more_complicated_inverse_linear_operator!(y,x)

# Arguments

- A: (array). Array used to define the linear operator

# Return

even_more_complicated_linear_operator!(y,x)

# Description

- Testing columnwise solve and reduction for an operator

- Acts on a vector "a" with size(a) = (n_i, n_j, n_k, n_s, n_ev, n_eh, n_r)

- This is a columnwise linear operator

# Arguments
- `y`: (array). OVERWRITTEN. output of linear operator acting on x
- `x`: (array). thing that the linear operator acts on
"""
function closure_even_more_complicated_linear_operator!(A)
  function even_more_complicated_linear_operator!(y,x)
    tup = size(y)
    for i6 in 1:tup[6]
      for i4 in 1:tup[4]
        for i2 in 1:tup[2]
          for i1 in 1:tup[1]
            tmp = x[i1, i2, :, i4, :, i6][:]
            tmp2 = A[i1,i2,i4,i6] * tmp
            y[i1, i2, :, i4, :, i6] .= reshape(tmp2, (tup[3], tup[5]))
          end
        end
      end
    end
  end
end

"""
closure_even_more_complicated_inverse_linear_operator!(A)

# Description
- Closure for even_more_complicated_inverse_linear_operator!(y,x)

# Arguments

- A: (array). Array used to define the linear operator

# Return

even_more_complicated_inverse_linear_operator!(y,x)
# Description

- Testing columnwise solve and reduction for an operator

- Acts on a vector "a" with size(a) = (n_i, n_j, n_k, n_s, n_ev, n_eh, n_r)

- This is a columnwise linear operator

# Arguments
- `y`: (array). OVERWRITTEN. output of linear operator acting on x
- `x`: (array). thing that the linear operator acts on
"""
function closure_even_more_complicated_inverse_linear_operator!(inv_A)
  function even_more_complicated_inverse_linear_operator!(y,x)
    tup = size(y)
    for i6 in 1:tup[6]
      for i4 in 1:tup[4]
        for i2 in 1:tup[2]
          for i1 in 1:tup[1]
            tmp = x[i1, i2, :, i4, :, i6][:]
            tmp2 = inv_A[i1,i2,i4,i6] * tmp
            y[i1, i2, :, i4, :, i6] .= reshape(tmp2, (tup[3], tup[5]))
          end
        end
      end
    end
  end
end


# Define the functions
tup = (3,4,5, 1, 2, 1)
y = randn(tup)
x = copy(y)
z = copy(x)

A = [randn(tup[3]*tup[5], tup[3]*tup[5]) for i1 in 1:tup[1], i2 in 1:tup[2], i4 in 1:tup[4], i6 in 1:tup[6] ]
# make sure it is positive definite
A = [A[i1,i2,i4,i6] * A[i1,i2,i4,i6]' + I for i1 in 1:tup[1], i2 in 1:tup[2], i4 in 1:tup[4], i6 in 1:tup[6] ]
inv_A = [inv(A[i1,i2,i4,i6]) for i1 in 1:tup[1], i2 in 1:tup[2], i4 in 1:tup[4], i6 in 1:tup[6] ]

new_linear_operator! = closure_even_more_complicated_linear_operator!(A)
new_inverse_linear_operator! = closure_even_more_complicated_inverse_linear_operator!(inv_A)


y = randn(tup)
x = copy(y)
z = copy(x)

new_linear_operator!(y,x)
new_inverse_linear_operator!(x, y)

norm(z .- x)

method(x,tol) = ColumnwisePreconditionedConjugateGradient(x, max_iter= tup[3]*tup[5] + 10, dims = (3,5))
linearsolver = method(x, tol)

y = copy(x)
iters = linearsolve!(new_linear_operator!, linearsolver, x, y; max_iters=Inf)
exact = copy(x)
new_inverse_linear_operator!(exact, y)
println("the norm of the error is ")
println(norm(x - exact)/norm(exact))
println("The number of iterations is ")
println(iters)
