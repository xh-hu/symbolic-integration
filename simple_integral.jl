# Testing file for integration
# Always use LiteralReal instead of Real!

module SymbolicIntegration

using Symbolics
using SymbolicUtils

export integral, integrate, integrate_from_table
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


# Attempt at writing an integration function
let
    INTEGRATION_TABLE = [
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

    global integrate_from_table
    global integrate

    function in_integration_table(op)
        return op == log || op == exp || op == sin || op == cos || op == tan || op == sinh || op == cosh || op == tanh || op == ^
    end

    function integrate_from_table(x)
        # @show x
        @assert (operation(x) == integral) "can only integrate an integral term"
        rule_tree = SymbolicUtils.Fixpoint(SymbolicUtils.Chain(INTEGRATION_TABLE))
        # Check if integration worked
        @assert (operation(rule_tree(x)) != integral) "integration expression is not currently supported"
        rule_tree(x)
    end

    """
    Integrates an expression x with only one variable
    """
    function integrate(x)
        # @show SymbolicUtils.inspect(x)
        # @show arguments(x)
        @assert (operation(x) == integral) "can only integrate an integral term"
        @assert (length(arguments(x)) == 2) "format: integral(x, y), which represents an integral of expression x with respect to y"
        e = simplify(arguments(x)[1], expand=true)
        v = simplify(arguments(x)[2])
        if is_exp_with_op(e, +)
            res = 0
            # @show arguments(e)
            for arg in arguments(e)
                # @show res
                res += integrate(integral(arg, v))
            end
            return res
        elseif is_exp_with_op(e, -) && size(arguments(e)) == 2
            return integrate(integral(arguments(e)[1], v)) - integrate(integral(arguments(e)[2], v))
        elseif is_exp_with_op(e, -) && size(arguments(e)) == 1
            return -integrate(integral(arguments(e)[1], v))
        elseif is_exp_with_op(e, *) && SymbolicUtils.is_literal_number(arguments(e)[1])
            # @show e
            return arguments(e)[1] * integrate(integral(arguments(e)[2], v))
        else
            for factor in factorize(e)
                # @show factor
                k = simplify_fractions(e/(factor * expand_derivatives(Differential(v)(factor))))
                # k = divides_to(e, factor * expand_derivatives(Differential(v)(factor)))
                # @show k
                if free_of(k, v)
                    if is_exp_with_op(factor, ^)
                        u = arguments(factor)[1]
                        n = arguments(factor)[2]
                        @assert SymbolicUtils._isinteger(n)
                        if n === -1
                            return k*log(u)
                        else
                            return k*u^(n+1)/(n+1)
                        end
                    else
                        return factor * v
                    end
                elseif istree(factor) && in_integration_table(operation(factor))
                    y = arguments(factor)[1]
                    # @show y
                    k2 = simplify_fractions(e/(factor * expand_derivatives(Differential(v)(y))))
                    # @show k2
                    if free_of(k2, v)
                        return k2*integrate_from_table(integral(factor, y))
                    end
                end
            end
            throw(error("integration expression is not currently supported"))
        end
    end

end

# @syms x::Real
# @show integrate(integral(cos(3*x)/(1-sin(3*x))^2, x))

end