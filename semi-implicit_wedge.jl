using GLMakie, XLSX
using LinearAlgebra


M = 2.0        
m = 1.0   
g = 9.8      
θ = 30.0      
θ_rad = deg2rad(θ)

M_mat = [M+m          m*cos(θ_rad);
         m*cos(θ_rad) m           ]

F_vec = [0.0;
         m * g * sin(θ_rad)]

a_exact = M_mat \ F_vec
a_X = a_exact[1]
a_s = a_exact[2]


T_end = 2.0
Δt = 0.0625
t = 0.0:Δt:T_end
nₜ = length(t)

q = zeros(2, nₜ)
v = zeros(2, nₜ)   


for n in 1:nₜ-1
    a_n = M_mat \ F_vec
    v[:, n+1] = v[:, n] + Δt * a_n
    q[:, n+1] = q[:, n] + Δt * v[:, n]
end


q_exact = zeros(2, nₜ)
for n in 1:nₜ
    q_exact[1, n] = 0.5 * a_X * t[n]^2
    q_exact[2, n] = 0.5 * a_s * t[n]^2
end

e = q - q_exact
# e_X_abs = abs.(e[1, :])
# e_s_abs = abs.(e[2, :])
e_X_err = e[1, :]
e_s_err = e[2, :]

fig = Figure(size = (1000, 800))


ax1 = Axis(fig[1, 1],
    title = "Displacement of 2-DOF Wedge-Block System (Semi-implicit Euler)",
    xlabel = "Time (s)", ylabel = "Displacement (m)")

lines!(ax1, t, q[1, :], color = :blue, linewidth = 2, label = "X_num (Wedge M)")
lines!(ax1, t, q[2, :], color = :red, linewidth = 2, label = "s_num (Block m rel-disp)")
lines!(ax1, t, q_exact[1, :], color = :gray, linestyle = :dash, linewidth = 2, label = "Exact Solution (X&s)")
lines!(ax1, t, q_exact[2, :], color = :gray, linestyle = :dash, linewidth = 2)

axislegend(ax1, position = :lt)
xlims!(ax1, 0, T_end)

ax2 = Axis(fig[2, 1],
    title = "Displacement Error Magnitude (|num - exact|)",
    xlabel = "Time (s)", ylabel = "Error |e| (m)")

lines!(ax2, t, e_X_err, color = :blue, linewidth = 2, label = "|X_num - X_exact|")
lines!(ax2, t, e_s_err, color = :red, linewidth = 2, label = "|s_num - s_exact|")

axislegend(ax2, position = :lt)
xlims!(ax2, 0, T_end)

fig

# save("d:/hmd/fig/一维/wedge_block_2dof.png", fig)

# XLSX.openxlsx("./excel/双自由度.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     for i in 1:length(t)
#         Sheet["Y$(i)"] = t[i]
#         Sheet["Z$(i)"] = q[1, i]
#         Sheet["AA$(i)"] = q[2, i]
#         Sheet["AB$(i)"] = q_exact[1, i]
#         Sheet["AC$(i)"] = q_exact[2, i]
#         Sheet["AD$(i)"] = e[1, i]
#         Sheet["AE$(i)"] = e[2, i]
#     end
# end
