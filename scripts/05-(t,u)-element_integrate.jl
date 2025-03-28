using DrWatson
@quickactivate "masters2semester_project"

using Plots
using DynamicalSystems
using OrdinaryDiffEq: Tsit5

include(srcdir("NaraharaPapers.jl"))
using .NaraharaPapers

include(srcdir("MyIntegration.jl"))
using .MyIntegration

include(srcdir("config.jl"))
include(srcdir("tools.jl"))


#########################################################################################
#########################################################################################

# Set dimentional system parameters
R, L, C, V_B = 5, 5e-6, 2e-9, 0.2 # Ohm, H, F, V
R₀ = NaraharaPapers.V_p/NaraharaPapers.I_p
E = V_B/NaraharaPapers.I_p/R
δ = R₀/R
ε = C*R₀*R/L
p = (ε, E, δ)

# Set initial state
u₀, v₀ = 0.0, 0.0
U₀ = [u₀, v₀]

# Set integration span
t_start, t_end = 0.0, 10.0
t_span = (t_start, t_end)

# Import I-V curve data
IV_curve_u_data = NaraharaPapers.IV_curve_u_data
IV_curve_v_data = NaraharaPapers.IV_curve_v_data


#########################################################################################
#########################################################################################


# Choose approximation type, from "data_linear" or "data_poly"
approx_type = "data_linear"

# Create f(u) based in the approx_type
if approx_type == "data_linear"
    IV_curve_func(u) = linear_approx_symmetrical(u, IV_curve_u_data, IV_curve_v_data)
    
    # Set title for plot
    curr_title = "Integration, linear interpolation approx."
elseif approx_type == "data_poly"
    # Set polynomial coefficients for dimentional I-V curve poly approx
    a₀ = -0.0014491777067221945
    a₁ = 0.4047903610420292
    a₂ = -8.814345521066247
    a₃ = 94.05841173048864
    a₄ = -593.8071165028581
    a₅ = 2362.6810193881056
    a₆ = -5964.114570853057
    a₇ = 9209.238555278202
    a₈ = -7889.603353342613
    a₉ = 2861.212170878973
    as = [a₀,a₁,a₂,a₃,a₄,a₅,a₆,a₇,a₈,a₉]

    # Convert dimentional coefficients to dimentionalless
    bs = NaraharaPapers.convert_dim_to_dimless_coeff(as, NaraharaPapers.V_p, NaraharaPapers.I_p)

    # Create a poly approx function
    IV_curve_func(u) = polynom_approx(u, bs)

    # Set title for plot
    curr_title = "Integration, polynomial approx."
else
    error("Unknown `approx_type`")
end


#########################################################################################
#########################################################################################


model! = MyIntegration.create_element_dimless_function(IV_curve_func)
diffeq = (alg = Tsit5(), abstol = 1e-5, reltol = 1e-5)
system = CoupledODEs(model!, U₀, p; diffeq)
Y, t = trajectory(system, t_end; Δt=t_end/10000)

display(Y[1,:])

#########################################################################################
#########################################################################################


p1 = plot(title=curr_title,
    titlefontsize=10,
    xlabel="t", ylabel="u", 
    label=nothing,
    xlims=(minimum(t), maximum(t)),
    ylims=(minimum(Y[:,1]), maximum(Y[:,1])),
    minorgrid=true
)

# Coords
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# Integrated trajectory
plot!(p1, t, Y[:,1], ls=:solid, label="")

p = plot(p1, legend=false, size=(800,500), dpi=300, left_margin=(20,:px))

savedir = plotsdir(GENERAL_PLOT_SAVE_DIR, "05-(t,u)-element_integrate")
savename = "05-$(approx_type)-(t,u)-element_integrate_$(time_ns()).pdf"
savepath = joinpath(savedir, savename)
mkpath(savedir)
savefig(p, savepath)