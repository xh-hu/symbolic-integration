using Symbolics
using SymbolicUtils
using SymbolicUtils.Code

integral(x...) = term(integral, Symbolics.unwrap.(x)...)

# Returns true if e is an expression with operation op
function is_exp_with_op(e, op)
    return istree(e) && op === operation(e)
end

let
    STANDARD_FORMS = [
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

    global slagle_integrate

    @enum EdgeType begin
        AND = 1
        OR = 2
        LEAF = 3
    end

    mutable struct GoalNode
        expr::Any
        dvar::Any
        result::Union{Any, Nothing}
        # stores the function to compute its own reuslt from child results
        compute::Union{Function, Nothing}
        relevant::Int

        parent::Union{GoalNode, Nothing}
        # childtype describes the AND/OR of branches from this node to its children
        childtype::EdgeType
        children::Vector{GoalNode}
    end

    """
    Helper function: finds the depth of an expression. Used to calculate "cost" of attempt. (tested to work)
    """
    function get_depth(e, depth)
        if !istree(e)
            return depth
        else
            res = 0
            for arg in arguments(e)
                res = max(res, get_depth(arg, depth+1))
            end
            return res
        end
    end

    function char_lt(x, y)
        return get_depth(x, 0) < get_depth(y, 0)
    end
    
    """
    Debugging function: prints the entire goal tree, by printing info from all nodes in prefix order
    """
    function get_root(node)
        if node.parent !== nothing
            return get_root(node.parent)
        end
        return node
    end

    function print_subtree(node)
        print(node)
        for child in node.children
            print_tree(child)
        end
    end

    # input is any node in the tree
    function print_tree(node)
        g = get_root(node)
        print_subtree(g)
    end

    """
    Helper function: prunes the goal tree. Takes in a GoalNode as input.
    """
    function prune_tree(node)
        println(node)
        node.relevant = 0
        # if the node is the root node (original goal), return true
        if node.parent === nothing
            return true
        end
        # prune all its children
        for child in node.children
            # For now, we do not check if node has multiple parents, since it is not implemented yet
            if child.relevance !== 0
                prune_tree(child)
            end
        end
        # prune ancestors if achievable
        if node.result !== nothing
            p = node.parent
            if p.childtype === OR::EdgeType
                # COMPUTE RESULT FROM CHILDREN
                p.result = p.compute(node.result)
                @show node.result
                @show p.result
                prune_tree(p)
            elseif p.childtype === AND::EdgeType
                total_rel = 0
                children_results = []
                for child in p.children
                    total_rel += child.relevance
                    push!(children_results, child.result)
                end
                if total_rel === 0
                    p.result = p.compute(children_results...)
                    prune_tree(p)
                end
            end
        end
    end

    """
    imsln: applies immediate solving to a goal list and mutates temp_goal_list input to be the new temporary goal list
    """
    function imsln(goal_list, temp_goal_list)
        if length(goal_list) === 0
            return false
        else
            g1 = goal_list[1]
            println(g1)
            if g1.relevant === 0
                goal_list = goal_list[2:end]
                imsln(goal_list)
            else
                # we skip c: reassigning parents for same expressions for now for simplicity
                # check to achieve g1: use simple integral temporarily
                rule_tree = SymbolicUtils.Fixpoint(SymbolicUtils.Chain(STANDARD_FORMS))
                if operation(rule_tree(integral(g1.expr, g1.dvar))) != integral
                    g1.result = rule_tree(integral(g1.expr, g1.dvar))
                    if prune_tree(g1) == true
                        return true
                    else
                        goal_list = goal_list[2:end]
                        imsln(goal_list, temp_goal_list)
                    end
                else
                    # Algorithmic transformations
                    if is_exp_with_op(g1.expr, *) && SymbolicUtils.is_literal_number(arguments(g1.expr)[1])
                        function f(x)
                            return  arguments(g1.expr)[1] * x
                        end
                        g1.compute = f
                        g2 = GoalNode(arguments(g1.expr)[2], g1.dvar, nothing, nothing, 1, g1, LEAF, [])
                        g1.childtype = OR
                        push!(g1.children, g2)
                        push!(goal_list, g2)
                    else
                        push!(temp_goal_list, g1)
                    end
                    goal_list = goal_list[2:end]
                    imsln(goal_list, temp_goal_list)
                end
            end
        end

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

        original_goal = GoalNode(e, v, nothing, nothing, 1, nothing, LEAF, [])
        goal_list = [original_goal];
        temp_goal_list = [];
        heuristic_list = [];

        if imsln(goal_list, temp_goal_list)
            return original_goal.result
        end

        # check resource allotment (maybe some timed condition?)
        while true
            heuristic_list = sort(temp_goal_list, lt=char_lt)
            temp_goal_list = []

            if length(heuristic_list) === 0
                throw(error("integration expression is not currently supported"))
            end

            g_i = heuristic_list[1]
            # for each heuristic transformation
            # if applicable
            g = g_i # apply to g_i
            g_i.childtype = LEAF
            push!(g_i.children, g)
            push!(goal_list, g)
            if imsln(goal_list, temp_goal_list)
                return original_goal.result
            else
                if g_i.result != nothing
                    continue
                else
                    break
                end
            end
            # end for loop
        end
    end

end

@syms x::Complex
@show slagle_integrate(integral(5*x^3, x))