include("./simple_integral.jl")
using .SymbolicIntegration
using Symbolics
using SymbolicUtils

@syms x::Complex

# Testing integration function:
# One term, directly from integral table
@show integrate(integral(exp(x), x))
# 2 terms, directly from integral table, using + and constant multiplier
@show integrate(integral(2*sin(x) + cos(x), x))
# 3 terms, directly from integral table, using -
@show integrate(integral(tan(x) + term(log, x) - sinh(x), x))
# Using derviative divides method only
@show integrate(integral(x*cos(x^2), x))
# Using both derivative divides and integral table rules
@show integrate(integral(cos(x)tan(2*sin(x)) + 3*exp(x) - x^3, x))
# Integrating on not integral - should throw error
try
    @show integrate(x^2+2*x)
catch e
    println("Threw error successfully")
end
# Integral does not have two terms - should throw error
try
    @show integrate(integral(sin(x)))
catch e
    println("Threw error successfully")
end
# Not supported - should throw error
try
    @show integrate(integral(x^2*sin(2*x)))
catch e
    println("Threw error successfully")
end