using Plots
using Images
using Statistics
using DelimitedFiles

# 初期条件
windV = 0.0 #[m/s]
initialV = 6.38 # [m/s]
initialT = 0.0 # [s]
initialRe = 0.001 # Reynolds number
timestep = 0.001 # dt [s]

s_dp = 20.0
in_dp = 1.0
e_dp = 300.0
rho_p = 1000 #kg/m^3 粒子密度
rho_f = 1.2 #kg/m^3 流体密度
mu = 1.8 * 1e-5 #Pa•s = N / m^2 • s 流体粘性係数
g = 9.807 # kg m/sec^2
max_time = 10.0

ZEP = 1e-8
MIN_RE = 1e-8 # Reynolds数の最小値

# Runge-Kutta法の関数
function rungeKutta4(dp, v, t, dt, cr, f)
    dpm = dp
    k1 = dt * f(dpm, t, v, cr)
    k2 = dt * f(dpm, t + dt/2.0, v + k1/2.0, cr)
    k3 = dt * f(dpm, t + dt/2.0, v + k2/2.0, cr)
    k4 = dt * f(dpm, t + dt, v + k3, cr)
    return v + (k1 + 2*k2 + 2*k3 + k4) / 6.0
end

# ODEの関数
function ODE(dp, t, v, cr)
    a = g 
    b = (3.0 * rho_f) / (4.0 * rho_p * dp * 1e-6)
    # return a - b * cr * (v - windV) * (v - windV)
    return a - b * cr * v * v
end

function simulate(dp, initialV, initialT, initialRe, timestep, ZEP, max_time)
    v = initialV - windV # m/s
    t = initialT # s
    dpm = dp
    if v == 0
        re = 0.001
    else 
        re = abs(v) * dpm * rho_f / mu * 1e-6 
    end
    
    cr = (2.25 * (re^(-0.31)) + 0.36 * (re^(0.06)))^3.45
    
    distance = 0.0 # m
    results = []
    
    max_re = re
    min_re = re

    while t <= max_time
        if isnan(v) || isinf(v) || v < -1e6  # 発散チェック
            println("Warning: Velocity diverged for dp=$dp um at time $t, skipping this data point.")
            return [], 0.0, NaN, NaN  # 発散した場合は空データを返す
        end

        if v > 0
            distance += v * timestep
        end

        newV = rungeKutta4(dpm, v, t, timestep, cr, ODE)
        re = abs(newV) * dpm * rho_f / mu * 1e-6  
        
        max_re = max(max_re, re)
        min_re = min(min_re, re)
        
        cr = (2.25 * (re^(-0.31)) + 0.36 * (re^(0.06)))^3.45
        push!(results, (t, v, distance))  
        t += timestep

        if abs(newV - v) <= ZEP
            v = newV  
            break
        end

        v = newV
    end

    if t > max_time
        println("Simulation reached max_time without full convergence.")
    end

    println("Final Velocity = $v for dp=$dp um")
    println("Fall Distance = $(distance * 100) mm ($dp um)")
    println("Reynolds Number Range for dp=$dp um: Min = $min_re, Max = $max_re")
    return results, distance, max_re, min_re
end

function main()
    par_size = []
    f_vel = []
    all_max_re = []
    all_min_re = []

    for idx in s_dp:in_dp:e_dp
        results, _, max_re, min_re = simulate(idx, initialV, initialT, initialRe, timestep, ZEP, max_time)
        
        if !isempty(results) && !isnan(results[end][2])  # 発散したデータを無視
            final_velocity = results[end][2]  
            # push!(f_vel, (final_velocity + windV) * 10e5 / 4000)
            push!(f_vel, (final_velocity + windV))
            push!(par_size, idx)
            push!(all_max_re, max_re)
            push!(all_min_re, min_re)
        else
            println("Skipping dp=$idx um due to invalid final velocity.")
        end
    end

    output_data = hcat(par_size, f_vel)
    writedlm("terminal_velocity-noflowpow.dat", output_data, ',')
    default(fontfamily = "serif-roman")
    plot(par_size, f_vel,
        st=:line, 
        lw=2,
        color=:black,
        xlabel="Particle [μm]", 
        ylabel="Final Velocity [m/s]", 
        guidefontsize=16,   # 軸ラベルのフォントサイズ
        tickfontsize=12,    # 軸の目盛りフォントサイズ
        # title="Particle Size vs Final Velocity", 
        legend=false,
        # title=:none,
        framestyle=:box,
        grid=false,
        )
    savefig("terminalvel-noflowpow.png")
end

main()
