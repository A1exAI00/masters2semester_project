using DrWatson
@quickactivate "masters2semester_project"

using Plots

include(srcdir("NaraharaPapers.jl"))
using .NaraharaPapers

include(srcdir("config.jl"))
include(srcdir("tools.jl"))


#########################################################################################
#########################################################################################


# Detect a number of threads; used to parallelize calculation
N_threads = Threads.nthreads()
println("Threads: $(N_threads)")

# Import I-V curve data
target_u_data = NaraharaPapers.IV_curve_u_data
target_v_data = NaraharaPapers.IV_curve_v_data

# Boundries for u
u_start, u_end, u_N = minimum(target_u_data), 1.9*maximum(target_u_data), 200
u_range = range(u_start, u_end, u_N)

# Boundries for root finding
u_start_root, u_end_root = u_start, u_end
u_span_root = (u_start_root, u_end_root)
u_eps_root = 1e-2


#########################################################################################
#########################################################################################


# Choose approximation type, from "data_linear" or "data_poly"
approx_type = "data_poly"

# Create f(u) based in the approx_type
if approx_type == "data_linear"
    # Create linear interpolation approx function
    f(u) = linear_approx_symmetrical(u, target_u_data, target_v_data)

    # Set title for plot
    curr_title = "Number of states of equalibrium, linear interpolation approx."
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
    f(x) = polynom_approx(x, bs)

    # Set title for plot
    curr_title = "Number of states of equalibrium, polynomial approx."
else
    error("Unknown `approx_type`")
end


#########################################################################################
#########################################################################################


# Create other function (you know what i mean)
g(u, E, δ) = E - δ*u

# Boundries for parameter E
E_start, E_end, E_N = 0.1, 4.0, 500
E_range = range(E_start, E_end, E_N)

# Boundries for parameter δ
δ_start, δ_end, δ_N = 0.0, 3.0, 500
δ_range = range(δ_start, δ_end, δ_N)

# Create a function to find a number of roots for specific parameters
function N_eq_func(E, δ)
    return length(find_roots(f, x->g(x, E, δ), u_span_root, u_eps_root))
end

# Iterate over a grid of δ and E values
eq_N = zeros(δ_N, E_N)
Threads.@threads for i in eachindex(E_range)
    if i < E_N/N_threads print("$(i) / $(Int(round(E_N/N_threads)))" * " "^10 * "\r") end
    E = E_range[i]
    for j in eachindex(δ_range)
        δ = δ_range[j]
        eq_N[j,i] = N_eq_func(E, δ)
    end
end


#########################################################################################
#########################################################################################


p1 = plot(title=curr_title,
    titlefontsize=10,
    xlabel="E", ylabel="δ", 
    label=nothing,
    xlims=(E_start, E_end),
    ylims=(δ_start, δ_end),
    minorgrid=true
)

# Heatmap with number of SoE denoted by color
heatmap!(p1, E_range, δ_range, eq_N)

p = plot(p1, legend=false, size=(800,500), dpi=300, left_margin=(20,:px))

savedir = plotsdir(GENERAL_PLOT_SAVE_DIR, "01-(E,δ)_num_of_SoE")
savename = "01-$(approx_type)-(E,δ)_num_of_SoE_$(time_ns()).pdf"
savepath = joinpath(savedir, savename)
mkpath(savedir)
savefig(p, savepath)