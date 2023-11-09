include("./SymbolicIntegration.jl")
using .SymbolicIntegration
using Symbolics
using SymbolicUtils

@syms x::Complex

# Testing integration function:
# One term, directly from integral table
@show simple_integrate(integral(exp(x), x))
# 2 terms, directly from integral table, using + and constant multiplier
@show simple_integrate(integral(2*sin(x) + cos(x), x))
# 3 terms, directly from integral table, using -
@show simple_integrate(integral(tan(2*x) + term(log, x) - sinh(x), x))
# Using derviative divides method only
@show simple_integrate(integral(x*cos(x^2), x))
# Using both derivative divides and integral table rules
@show simple_integrate(integral(cos(x)tan(2*sin(x)) + 3*exp(x) - x^3, x))
# Integrating on not integral - should throw error
try
    @show simple_integrate(x^2+2*x)
catch e
    println("Threw error successfully")
end
# Integral does not have two terms - should throw error
try
    @show simple_integrate(integral(sin(x)))
catch e
    println("Threw error successfully")
end
# Not supported - should throw error
try
    @show simple_integrate(integral(x^2*sin(2*x)))
catch e
    println("Threw error successfully")
end