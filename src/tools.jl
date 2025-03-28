# For find_root
using NonlinearSolve: IntervalNonlinearProblem, solve

#########################################################################################
#########################################################################################

"""
    linear_approx(x, xs, ys, k_before, k_after)

Function to perform a linear interpolation for `x∈(min(xs), max(xs))` and linear exterpolation for `x∉(min(xs), max(xs))`

Parameters:

`x`: point to approximate

`xs`: tuple/array of known X coordinates (in order, no duplecates)

`ys`: tuple/array of known Y coordinates (in order, no duplecates)

`k_before`: slope before the first point (to extrapolate)

`k_after`: slope after the last point (to extrapolate)
"""
function linear_approx(x, xs, ys, k_before, k_after)
    # if need to exterpolate
    if x ≤ xs[1]
        return k_before*(x-xs[1]) + ys[1]
    end
    if x ≥ xs[end]
        return  k_after*(x-xs[end]) + ys[end]
    end

    # if need to interpolate
    for i in 1:length(xs)-1
        if xs[i] ≤ x ≤ xs[i+1]
            k_curr = (ys[i+1]-ys[i])/(xs[i+1]-xs[i])
            return k_curr*(x-xs[i]) + ys[i]
        end
    end
end

"""
    linear_approx_after(x, xs, ys, k)

Function to perform a linear interpolation for `x∈(min(xs), max(xs))` and linear exterpolation for `x∉(min(xs), max(xs))`

Parameters:

`x`: point to approximate

`xs`: tuple/array of known X coordinates (in order, no duplecates)

`ys`: tuple/array of known Y coordinates (in order, no duplecates)

`k`: slope after the last point (to extrapolate)
"""
function linear_approx_after(x, xs, ys, k)
    k_before = (ys[2]-ys[1])/(xs[2]-xs[1])
    return linear_approx(x, xs, ys, k_before, k)
end

"""
    linear_approx_symmetrical(x, xs, ys)

Function to perform a linear interpolation for `x∈(min(xs), max(xs))` and linear exterpolation for `x∉(min(xs), max(xs))`

Extrapolation is done by continuing the first and the last line segments.

Parameters:

`x`: point to approximate

`xs`: tuple/array of known X coordinates (in order, no duplecates)

`ys`: tuple/array of known Y coordinates (in order, no duplecates)
"""
function linear_approx_symmetrical(x, xs, ys)
    k_before = (ys[2]-ys[1])/(xs[2]-xs[1])
    end_i = length(xs)
    k_after = (ys[end_i]-ys[end_i-1])/(xs[end_i]-xs[end_i-1])
    return linear_approx(x, xs, ys, k_before, k_after)
end

#########################################################################################
#########################################################################################

"""
    polynom_approx(x, as)

Function to calculate a value of polynomial with coefficients `as` at `x`.

Parameters:

`x`: point to approximate

`as`: tuple/array of polynomial coefficients: (a₀, a₁, ...)
"""
function polynom_approx(x, as)
    res = 0.0
    for i in 0:length(as)-1
        res += as[i+1]*x^i
    end
    return res 
end

#########################################################################################
#########################################################################################

"""
    sub_index(i_char)

Function to convert integer charecter `i_char` to sub index

Non integer charecters are ignored.

Parameters:

`i_char`: integer charecter
"""
function sub_index(i_char)
    if i_char == '1' return "₁"
    elseif i_char == '2' return "₂"
    elseif i_char == '3' return "₃"
    elseif i_char == '4' return "₄"
    elseif i_char == '5' return "₅"
    elseif i_char == '6' return "₆"
    elseif i_char == '7' return "₇"
    elseif i_char == '8' return "₈"
    elseif i_char == '9' return "₉"
    elseif i_char == '0' return "₀"
    else return ""
    end
end

"""
    sub_string(index::String)

Function to convert a string of integer charecters to a string of sub indexes

Non integer charecters are ignored.

Parameters:

`index::String`: a string of integer charecters
"""
function sub_string(index::String)
    new_string = ""
    for i in eachindex(index)
        char = index[i]
        # println(char)
        # println(sub_index(char))
        new_string *= sub_index(char)
    end
    return new_string
end

"""
    sub_string(index::Int)

Function to convert an integer to a string of sub index.

Parameters:

`index::Int`: an integer to convert to sub index
"""
function sub_string(index::Int)
    return sub_string(string(index))
end

#########################################################################################
#########################################################################################

"""
    local_maximum(arr)

Function to find local maximum values in an array

Parameters:

`arr`: array if comparable elements
"""
function local_maximum(arr)
    vals = []
    for i in eachindex(arr) 
        if i == 1
            continue
        elseif i == length(arr)
            continue
        end
        if (arr[i-1] < arr[i]) && (arr[i] > arr[i+1])
            push!(vals, arr[i])
        end
    end
    return vals
end

#########################################################################################
#########################################################################################

"""
    find_root(f, g, x_span)

Function to find a root `f(x)=g(x)` at specific `x_span`.

Parameters:

`f`: first function

`g`: second function

`x_span`: x window to search root in
"""
function find_root(f, g, x_span)
    delta(x, p) = f(x) - g(x)
    if delta(x_span[1], 0) * delta(x_span[2], 0) > 0
        return nothing
    end
    prob = IntervalNonlinearProblem(delta, x_span)
    sol = solve(prob)
    return sol.u
end

"""
    find_roots(f, g, x_span, eps)

Function to find multiple roots `f(x)=g(x)` at specific `x_span`.

Parameters:

`f`: first function

`g`: second function

`x_span`: x window to search root in

`eps`: minimal distance between two roots
"""
function find_roots(f, g, x_span, eps)
    x_segm_N = Int(round((x_span[2] - x_span[1])/eps))

    x_roots = []
    for i in 0:x_segm_N-1
        bound_l = x_span[1] + i*eps
        bound_r = x_span[1] + (i+1)*eps
        root = find_root(f, g, (bound_l, bound_r))
        if !isnothing(root)
            push!(x_roots, root)
        end
    end
    return unique!(x_roots)
end

