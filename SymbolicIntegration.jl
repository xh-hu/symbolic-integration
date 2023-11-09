module SymbolicIntegration

using Symbolics
using SymbolicUtils

export integral, simple_integrate, integrate_from_table, slagle_integrate
integral(x...) = term(integral, Symbolics.unwrap.(x)...)


# Returns true if e is an expression with operation op
function is_exp_with_op(e, op)
    return istree(e) && op === operation(e)
end

# Factorizes a polynomial expression e: returns a list of expressions which are multiplied/divided to get e
function factorize(e)
    e = simplify(e)
    if !istree(e)
        return [e]
    elseif operation(e) == *
        return collect(Iterators.flatten([factorize(arguments(e)[i]) for i=1:length(arguments(e))]))
    elseif operation(e) == /
        return collect(Iterators.flatten([factorize(arguments(e)[1]), factorize(arguments(e)[2]) .^ -1]))
    else
        return [e]
    end
end

# Divides the numer expression by the denom expression
function divides_to(numer, denom)
    n_factors = factorize(simplify(numer))
    d_factors = factorize(simplify(denom))
    # @show n_factors
    # @show d_factors
    common_fact = [i for i in n_factors, j in d_factors if isequal(i, j)]
    # @show common_fact
    res_numer = 1
    for i in n_factors
        if length(findall(x->isequal(x, i), common_fact)) === 0
            res_numer = res_numer * i
        end
    end
    # @show res_numer
    res_denom = 1
    for i in d_factors
        if length(findall(x->isequal(x, i), common_fact)) === 0
            res_denom = res_denom * i
        end
    end
    # @show res_denom
    return res_numer/res_denom
end

# Returns true if e does not contain the symbol x
function free_of(e, x)
    if !istree(e)
        return e !== x
    else
        if size(arguments(e)) == 1
            return free_of(arguments(e)[1], x)
        else
            return free_of(arguments(e)[1], x) && free_of(arguments(e)[2], x)
        end
    end
end

global INTEGRATION_TABLE = [
    @rule(integral((~x)^(~y::(y -> SymbolicUtils._isinteger(y) && y !== -1)), ~x) => ((~x)^(~y+1))/(~y+1))
    @rule(integral(exp(~x), ~x) => exp(~x))
    @rule(integral((~y::(y -> SymbolicUtils.is_literal_number(y)))^(~x), ~x) => ((~y)^(~x))//log(~y))
    @rule(integral(log(~x), ~x) => (~x) * log(~x) - (~x))
    @rule(integral(log((~y::(y -> SymbolicUtils.is_literal_number(y))), ~x), ~x) => (~x) * log(~y, ~x) - (~x)//log(~y))
    @rule(integral(sin(~x), ~x) => -cos(~x))
    @rule(integral(cos(~x), ~x) => sin(~x))
    @rule(integral(tan(~x), ~x) => -log(cos(~x)))
    @rule(integral(cot(~x), ~x) => log(sin(~x)))
    @rule(integral(sec(~x), ~x) => log(sec(~x) + tan(~x)))
    @rule(integral(csc(~x), ~x) => log(csc(~x) - cot(~x)))
    @rule(integral(asin(~x), ~x) => (~x) * asin(~x) + sqrt(1-(~x)^2))
    @rule(integral(acos(~x), ~x) => (~x) * acos(~x) - sqrt(1-(~x)^2))
    @rule(integral(atan(~x), ~x) => (~x) * atan(~x) - 0.5 * log(1+(~x)^2))
    @rule(integral(acot(~x), ~x) => (~x) * acot(~x) + 0.5 * log(1+(~x)^2))
    @rule(integral(asec(~x), ~x) => (~x) * asec(~x) - log(~x + sqrt((~x)^2 - 1)))
    @rule(integral(acsc(~x), ~x) => (~x) * acsc(~x) + log(~x + sqrt((~x)^2 - 1)))
    @rule(integral(sec(~x)^2, ~x) => tan(~x))
    @rule(integral(csc(~x)^2, ~x) => -cot(~x))
    @rule(integral(sinh(~x), ~x) => cosh(~x))
    @rule(integral(cosh(~x), ~x) => sinh(~x))
    @rule(integral(tanh(~x), ~x) => log(cosh(~x)))
]

include("./simple_integral.jl")
include("./slagle_integral.jl")

end