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
l, r, N = 5e-6, 5, 3
R₀ = NaraharaPapers.V_p/NaraharaPapers.I_p
E = V_B/NaraharaPapers.I_p/R
δ = R₀/R
ε = C*R₀*R/L
η = l*R/L/R₀
κ = r/R₀

# Set initial state
U₀ = zeros(3N) .+ (rand(3N).-0.5).*0.002

# Import I-V curve data
IV_curve_V_data = NaraharaPapers.IV_curve_V_data
IV_curve_I_data = NaraharaPapers.IV_curve_I_data
IV_curve_u_data = NaraharaPapers.IV_curve_u_data
IV_curve_v_data = NaraharaPapers.IV_curve_v_data


#########################################################################################
#########################################################################################


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

# Choose approximation type, from "data_linear" or "data_poly"
approx_type = "data_poly"

# Choose dimentional or dimentionalless model
dimentional = false

if dimentional && approx_type == "data_linear"
    IV_curve_func(V) = linear_approx_symmetrical(V, IV_curve_V_data, IV_curve_I_data)
    ode! = MyIntegration.create_loop_dim_ode(IV_curve_func)

elseif dimentional && approx_type == "data_poly"
    IV_curve_func(V) = polynom_approx(V, as)
    ode! = MyIntegration.create_loop_dim_ode(IV_curve_func)

elseif !dimentional && approx_type == "data_linear"
    IV_curve_func(u) = linear_approx_symmetrical(u, IV_curve_u_data, IV_curve_v_data)
    ode! = MyIntegration.create_loop_dimless_ode(IV_curve_func)

elseif !dimentional && approx_type == "data_poly"
    bs = NaraharaPapers.convert_dim_to_dimless_coeff(as, NaraharaPapers.V_p, NaraharaPapers.I_p)
    IV_curve_func(u) = polynom_approx(u, bs)
    ode! = MyIntegration.create_loop_dimless_ode(IV_curve_func)
else
    error("Unknown `approx_type`")
end

# Set xlabel and ylabel for plot
plot_xlabel = if dimentional "t" else "τ" end
plot_ylabel = if dimentional "Iₙ" else "vₙ" end

# Set a title for plot
plot_title = ""

# Diffetenr end times for different models
t_end = if dimentional 1e-5 else 10.0 end

p = if dimentional
    (R, L, C, V_B, l, r, N)
else 
    (ε, E, δ, η, κ, N)
end


#########################################################################################
#########################################################################################


diffeq = (alg = Tsit5(), abstol = 1e-5, reltol = 1e-5)
system = CoupledODEs(ode!, U₀, p; diffeq)
Y, t = trajectory(system, t_end; Δt=t_end/10000)


#########################################################################################
#########################################################################################


p1 = plot(title=curr_title,
    titlefontsize=10,
    xlabel=plot_xlabel, ylabel=plot_ylabel, 
    label=nothing,
    xlims=(minimum(t), maximum(t)),
    ylims=(minimum(Y[:,1]), maximum(Y[:,1])),
    minorgrid=true
)

# Coords
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# Plot for all oscillator element
for i in 1:N
    plot!(p1, t, Y[:,i], ls=:solid, label="")
end

p = plot(p1, legend=false, size=(800,500), dpi=300, left_margin=(20,:px))

savedir = plotsdir(GENERAL_PLOT_SAVE_DIR, "05-Oscill-loop_integrate")
savename = "05-$(approx_type)-$(dimentional)-Oscill-loop_integrate_$(time_ns()).pdf"
savepath = joinpath(savedir, savename)
mkpath(savedir)
savefig(p, savepath)