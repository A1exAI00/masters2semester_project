using DrWatson
@quickactivate "masters2semester_project"

using Plots
using Optim

include(srcdir("NaraharaPapers.jl"))
using .NaraharaPapers

include(srcdir("config.jl"))
include(srcdir("tools.jl"))


#########################################################################################
#########################################################################################


# Nonuniform I-V curve sampling
target_Vs_nonuniform = NaraharaPapers.IV_curve_V_data
target_Is_nonuniform = NaraharaPapers.IV_curve_I_data
curr_linear_approx(x) = linear_approx_symmetrical(x, target_Vs_nonuniform, target_Is_nonuniform)

# Interpolate I-V curve uniformly
N_uniform = 200
target_Vs_uniform = range(minimum(target_Vs_nonuniform), 1.1*maximum(target_Vs_nonuniform), N_uniform)
target_Is_uniform = curr_linear_approx.(target_Vs_uniform)

# Assign weights to points on I-V curve
high_weight_values = [20.0, 20.0, 7.0]
high_weight_ares = [
    (0.05, 0.065),
    (0.15, 0.34),
    (0.4, 0.45)
]

function calc_weights!(weights, high_weight_values, high_weight_ares)
    step_uniform = (maximum(target_Vs_nonuniform) - minimum(target_Vs_nonuniform))/N_uniform
    for i in eachindex(weights)
        for j in eachindex(high_weight_ares)
            area_ = high_weight_ares[j]
            w_ = high_weight_values[j]
            if area_[1] ≤ i*step_uniform ≤ area_[2]
                weights[i] = w_
            end
        end
    end
end

weights = ones(length(target_Vs_uniform))
calc_weights!(weights, high_weight_values, high_weight_ares)


#########################################################################################
#########################################################################################

# Optim parameters
OPTIM_G_TOL = 0.0
OPTIM_ITERATIONS = 100000
OPTIM_ALG = LBFGS()
OPTIM_DISPLAY = false


function general_poly_weighted_loss(params, x_data, y_data, weights)
    func_(x) = polynom_approx(x, params)
    sq_err = weights .* (y_data .- func_.(x_data)).^2
    return sum(sq_err)/length(sq_err)
end

function poly_optimize(params₀, x_data, y_data, weights)
    loss_(params) = general_poly_weighted_loss(params, x_data, y_data, weights)
    options = Optim.Options(g_tol=OPTIM_G_TOL, iterations=OPTIM_ITERATIONS)
    optim_result = optimize(loss_, params₀, OPTIM_ALG, options)
    if OPTIM_DISPLAY display(optim_result) end
    return Optim.minimizer(optim_result)
end

function easy_poly_optimize(degree, x_data, y_data, weights)
    return 
end


#########################################################################################
#########################################################################################

degree = 9

# Calc oprimized parameters and optimized poly function
optimized_params = poly_optimize(zeros(degree+1), target_Vs_uniform, target_Is_uniform, weights)
curr_poly_approx(x) = polynom_approx(x, optimized_params)

# Output coefficients
println("degree: $(degree)")
for i in eachindex(optimized_params)
    println("   a$(sub_string(string(i-1))) = $(optimized_params[i])")
end


#########################################################################################
#########################################################################################


p1 = plot(title="NEC1S1763 I-V curve polynomial approx.",
    titlefontsize=10,
    xlabel="V, В", ylabel="I, A", 
    label=nothing,
    xlims=(minimum(target_Vs_uniform), maximum(target_Vs_uniform)),
    ylims=(minimum(target_Is_uniform), maximum(target_Is_uniform)),
    minorgrid=true
)

# Coords
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# Experimental data
plot!(p1, target_Vs_uniform, target_Is_uniform, ls=:solid, label="Experimental data")

# Poly approx
plot!(p1, target_Vs_uniform, curr_poly_approx.(target_Vs_uniform), ls=:solid, label="Polynomial approx.")

p = plot(p1, legend=true, size=(800,500), dpi=300, left_margin=(20,:px))

savedir = plotsdir(GENERAL_PLOT_SAVE_DIR, "00-(V,I)-calc_polynom_as")
savename = "00-$(lpad(degree,3,"0"))-calc_polynom_as$(time_ns()).pdf"
savepath = joinpath(savedir, savename)
mkpath(savedir)
savefig(p, savepath)