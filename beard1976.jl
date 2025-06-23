using Printf
using DelimitedFiles

function terminal_vel(dp, visco, press, Kelvin)

    T = Kelvin #K
    p = press #mb
    η = visco * 0.1 #kg/m/s

    g = 9.817 #m/s2
    d_0 = dp * 1e-6 #m
    l_0 = 6.62e-6 * 0.01 #cm → m
    p_0 = 1013.25 #mb
    η_0 = 0.0001818 * 0.1 #g/cm/s → kg/m/s
    T_0 = 293.15 #K
    ρ_w = 1000 #kg/m3
    ρ_f = 1.225 #kg/m3
    Δρ = ρ_w - ρ_f #kg/m3

    if 0.5 <= dp && dp < 19.0
        C_l = Δρ * g / (18.0 * η) #1/m/s
        l = l_0 * (η / η_0) * (p_0 / p) * (T / T_0)^0.5 #m
        C_ac = 1 + 2.51 * l / d_0 #[-]
        return C_l * C_ac * d_0^2 #m/s
    elseif 19.0 <= dp && dp < 1070.0
        l = l_0 * (η / η_0) * (p_0 / p) * (T / T_0)^0.5 #m
        C_ac = 1 + 2.51 * l / d_0 #[-]
        C_2 = 4 * ρ_f * Δρ * g / (3.0 * η^2) #1/m3
        NDa = C_2 * d_0^3 #[-]
        X = log(NDa) #[-]
        Y = -0.318657e1 + 0.992696 * X -0.153193e-2 * X^2 -0.987059e-3 * X^3 -0.578878e-3 * X^4 + 0.855176e-4 * X^5 -0.327815e-5 * X^6 #[-]
        NRe = C_ac * exp(Y) #[-]
        return η * NRe / (ρ_f * d_0) #[m/s]
    elseif 1070.0 <= dp && dp < 7000.0
        println("Not implemented")
        return nothing
    else
        println("Outside the scope of the definition")
        return nothing
    end
end

#条件変更
η = 0.0001818 #g/cm/s
p = 1013.25 #mb
T = 293.15 #K

par_size =[]
ter_vel=[]
for idx in 20.0:1.0:300.0
    vel = terminal_vel(idx, η, p, T)
    push!(par_size, idx)
    push!(ter_vel, vel * 1e6 / 4000)
end

output_data = hcat(par_size, ter_vel)
writedlm("beard_tervel-300.dat", output_data, ',')
# dp = 30 #um

# #test
# dp_values = [10, 30, 100, 500, 1100, 8000]
# for idx in dp_values
#     println(terminal_vel(idx, η, p, T))
# end
