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
