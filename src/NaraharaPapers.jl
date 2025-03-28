module NaraharaPapers

using DrWatson
@quickactivate "masters2semester_project"

include(srcdir("tools.jl"))

#########################################################################################
#########################################################################################

function convert_V_to_u(V, V_p)
    return V/V_p
end

function convert_u_to_V(u, V_p)
    return u*V_p
end

function convert_I_to_v(I, I_p)
    return I/I_p
end

function convert_v_to_I(v, I_p)
    return v*I_p
end

function convert_dim_to_dimless_coeff(as, V_p, I_p)
    return [as[i]*V_p^(i-1)/I_p for i in eachindex(as)]
end

#########################################################################################
#########################################################################################

# Element parameters
# paper: Trigger waves in an oscillator loop and their applications in voltage-controlled oscillation with high tuning gain
L_se = L_sh = 5e-6 # H
C = 2e-9 # F
R_se = R_sh = 5 # Ohm

# I-V curve data for NEC 1S1763
# paper: Self-injection Locking of Rotary Traveling Pulses in Resonant-Tunneling-Diode Transmission-Line Loop
IV_curve_V_data = [0.0, 0.01028524, 0.02062652, 0.02902881, 0.03937009, 0.04971137, 0.06005265, 0.07039393, 0.07879622, 0.0891375, 0.09947878, 0.12016134, 0.13955124, 0.1602338, 0.17897737, 0.19965993, 0.24942734, 0.29984108, 0.34960849, 0.40002223, 0.44978964, 0.49955705, 0.54932446]
IV_curve_I_data = [0.0, 0.0016478873, 0.00299999994, 0.00396126752, 0.0047112675, 0.00519718298, 0.00548239424, 0.00557746466, 0.00353873232, 0.00329577458, 0.00320070416, 0.0030739436, 0.0030739436, 0.00292605628, 0.00288380276, 0.002862676, 0.00278873234, 0.00261971826, 0.0022288732, 0.00072887324, 0.00102464788, 0.00297887318, 0.00493309848]

#########################################################################################
#########################################################################################

I_p = local_maximum(IV_curve_I_data)[1]
I_p_index = findfirst(I->I==I_p, IV_curve_I_data)
V_p = IV_curve_V_data[I_p_index]

IV_curve_u_data = convert_V_to_u.(IV_curve_V_data, V_p)
IV_curve_v_data = convert_I_to_v.(IV_curve_I_data, I_p)

# I-V curve from data
I_D_func_data_dim(V) = linear_approx_symmetrical(V, IV_curve_V_data, IV_curve_I_data)
I_D_func_data_dimless(u) = linear_approx_symmetrical(u, IV_curve_u_data, IV_curve_v_data)

#########################################################################################
#########################################################################################

# function narahara_element_dim_data(dU, U, p, t)
#     L_sh, C, R_sh, V_B = p
#     V, I = u
#     du[1] = 1/C * (I - I_D_func_data(V))
#     du[2] = 1/L_sh * (V_B - I*R_sh - V)
# end

end