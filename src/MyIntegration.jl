module MyIntegration

using DynamicalSystems


#########################################################################################
#########################################################################################


"""
    create_element_dim_ode(BAX_dim::Function)

Function to create a model function for a dynamical system of one dimentional element with arbitrury I-V curve function.

Parameters:

`BAX_dim::Function`: an I-V curve function in dimentional variables
"""
function create_element_dim_ode(BAX_dim::Function)
    ode!(dU, U, p, t) = begin
        R, L, C, V_B = p
        V, I = U
        dU[1] = 1/C * (I - BAX_dim(V))
        dU[2] = 1/L * (V_B - I*R - V)
        return nothing
    end
    return ode!
end

"""
    create_element_dimless_ode(BAX_dimless::Function)

Function to create a model function for a dynamical system of one dimentionalles element with arbitrury I-V curve function.

Parameters:

`BAX_dimless::Function`: an I-V curve function in dimentionalless variables
"""
function create_element_dimless_ode(BAX_dimless::Function)
    ode!(dU, U, p, t) = begin
        ε, E, δ = p
        u, v = U
        dU[1] = 1/ε * (v - BAX_dimless(u))
        dU[2] = E - v - δ*u
        return nothing
    end
    return ode!
end


"""
    create_loop_dim_ode(BAX_dim::Function)

Function to create a model function for a dynamical system of a loop of coupled dimentional elements with arbitrury I-V curve function.

Parameters:

`BAX_dim::Function`: an I-V curve function in dimentional variables
"""
function create_loop_dim_ode(BAX_dim::Function)
    ode!(dU, U, p, t) = begin
        R, L, C, V_B, l, r, N = p
        Vₙ, Iₙ, iₙ = U[0N+1:1N], U[1N+1:2N], U[2N+1:3N]
        for i in 1:N
            prev_i = if i == 1 N else i-1 end
            next_i = if i == N 1 else i+1 end
            dU[0N+i] = 1/C * (Iₙ[i] - iₙ[i] + iₙ[prev_i] - BAX_dim(Vₙ[i]))
            dU[1N+i] = 1/L * (V_B - R*Iₙ[i] - Vₙ[i])
            dU[2N+i] = 1/l * (Vₙ[i] - Vₙ[next_i] - r*iₙ[i])
        end
        return nothing
    end
    return ode!
end


"""
    create_loop_dimless_ode(BAX_dim::Function)

Function to create a model function for a dynamical system of a loop of coupled dimentionalless elements with arbitrury I-V curve function.

Parameters:

`BAX_dim::Function`: an I-V curve function in dimentional variables
"""
function create_loop_dimless_ode(BAX_dimless::Function)
    ode!(dU, U, p, t) = begin
        ε, E, δ, η, κ, N = p
        uₙ, vₙ, wₙ = U[0N+1:1N], U[1N+1:2N], U[2N+1:3N]
        for i in 1:N
            prev_i = if i == 1 N else i-1 end
            next_i = if i == N 1 else i+1 end
            dU[0N+i] = 1/ε * (vₙ[i] - wₙ[i] + wₙ[prev_i] - BAX_dimless(uₙ[i]))
            dU[1N+i] = E - vₙ[i] - δ*uₙ[i]
            dU[2N+i] = 1/η * (uₙ[i] - uₙ[next_i] - κ*wₙ[i])
        end
        return nothing
    end
    return ode!
end


"""
    create_loop_resistive_dimless_ode(BAX_dim::Function)

Function to create a model function for a dynamical system of a loop of coupled dimentionalless elements with arbitrury I-V curve function.

Parameters:

`BAX_dim::Function`: an I-V curve function in dimentional variables
"""
function create_loop_resistive_dimless_ode(BAX_dimless::Function)
    ode!(dU, U, p, t) = begin
        ε, E, δ, κ, N = p
        uₙ, vₙ = U[0N+1:1N], U[1N+1:2N]
        for i in 1:N
            prev_i = if i == 1 N else i-1 end
            next_i = if i == N 1 else i+1 end
            dU[0N+i] = 1/ε * (vₙ[i] - BAX_dimless(uₙ[i]) - 1/κ * (uₙ[]))
            dU[1N+i] = E - vₙ[i] - δ*uₙ[i]
        end
        return nothing
    end
    return ode!
end

#########################################################################################
#########################################################################################

mutable struct Model
    ode!::Function
    variables_symb
    parameters_symb
    t_characteristic::Real
end


function integrate(model::Model, U₀, t_end, Δt)
    diffeq = (alg = Tsit5(), abstol = 1e-5, reltol = 1e-5)
    system = CoupledODEs(model.ode!, U₀, p; diffeq)
    Y, t = trajectory(system, t_end; Δt=Δt)
    return Y, t
end

#########################################################################################
#########################################################################################


# """
#     integrate(model, U₀, t_span, p; abstol=nothing, reltol=nothing, alg=nothing, save_at=nothing)

# Function to integrate model on a short time interval.

# Parameters:

# `model`: function of a model o integrate

# `U₀`: vector of initial state veriables

# `t_span`: tuple of start and end time stamps

# `p`: tuple of parameters

# `abstol`: absolute tolerance of integration

# `reltol`: reletive tolerance of integration

# `alg`: algorithm to perform an integration

# `save_at`: a vector of time stamps to save integration progress at
# """
# function integrate(model, U₀, t_span, p; abstol=nothing, reltol=nothing, alg=nothing, save_at=nothing)
#     ABSTOL = 1e-5
#     RELTOL = 1e-5
#     ALG = Tsit5()

#     if isnothing(reltol) reltol = RELTOL end
#     if isnothing(abstol) abstol = ABSTOL end
#     if isnothing(alg) alg = ALG end

#     prob = ODEProblem(model, U₀, t_span, p)
#     sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)
#     return sol
# end

# function integrate(model, U₀, t_span, p; diffeq=Nothing)
#     if isnothing(diffeq) diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9) end
#     system = ContinuousDynamicalSystem(model, U₀, p; diffeq)
#     return trajectory(system, t_span; Ttr = 2.2, Δt = sampling_time)
# end

end