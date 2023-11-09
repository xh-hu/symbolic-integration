# Attempt at writing slagle integration
let
    # global integrate_from_table
    global slagle_integrate

    function in_integration_table(op)
        return op == log || op == exp || op == sin || op == cos || op == tan || op == sinh || op == cosh || op == tanh || op == ^
    end

    """
    Integrates an expression x with only one variable
    """
    function slagle_integrate(x)
        # @show SymbolicUtils.inspect(x)
        # @show arguments(x)
        @assert (operation(x) == integral) "can only integrate an integral term"
        @assert (length(arguments(x)) == 2) "format: integral(x, y), which represents an integral of expression x with respect to y"
        e = simplify(arguments(x)[1], expand=true)
        v = simplify(arguments(x)[2])
        try
            return simple_integrate(x)
        catch e
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
                # goal tree impl here
                throw(error("integration expression is not currently supported"))
            end
        end
    end

end