# If using LiteralReal, change simplify_fractions in code to divides_to

# Attempt at writing an integration function
let
    global integrate_from_table
    global simple_integrate

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
    function simple_integrate(x)
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
                res += simple_integrate(integral(arg, v))
            end
            return res
        elseif is_exp_with_op(e, -) && size(arguments(e)) == 2
            return simple_integrate(integral(arguments(e)[1], v)) - simple_integrate(integral(arguments(e)[2], v))
        elseif is_exp_with_op(e, -) && size(arguments(e)) == 1
            return -simple_integrate(integral(arguments(e)[1], v))
        elseif is_exp_with_op(e, *) && SymbolicUtils.is_literal_number(arguments(e)[1])
            # @show e
            return arguments(e)[1] * simple_integrate(integral(arguments(e)[2], v))
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