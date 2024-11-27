using ApproxOperator, JuMP, Ipopt, CairoMakie, GLMakie, XLSX, LinearAlgebra

m = 1.0  
k = 1.0  
ω = sqrt(k/m)  
h = 0.01 
t_final = 10.0  

q₀ = 1.0  
p₀ = 0.0  

function evolve_one_step(q, p, h, m, k)
    q_new = q + h * p / m
    p_new = p - h * k * q
    return q_new, p_new
end

t = 0
q = q₀
p = p₀
q_history = [q]
p_history = [p]
time_history = [t]

while t < t_final
    q, p = evolve_one_step(q, p, h, m, k)
    push!(q_history, q)
    push!(p_history, p)
    push!(time_history, t)
    t += h
end

fig = Figure(resolution = (600, 400))

ax1 = Axis(fig, title = "Displacement vs Time")
lines!(ax1, time_history, q_history, label = "Displacement q(t)")

ax2 = Axis(fig, title = "Momentum vs Time")
lines!(ax2, time_history, p_history, label = "Momentum p(t)")

ax3 = Axis(fig, title = "Phase Space Trajectory")
scatter!(ax3, q_history, p_history, label = "Phase Space")

display(fig)

# face = zeros(nₑ,3)
# for (i,elm) in enumerate(elements["Ω"])
#     face[i,:] .= [x.𝐼 for x in elm.𝓒]
# end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
# meshscatter!(ax,xs,ys,ds,color=ds,markersize = 0.05)
# meshscatter!(ax,xs,ys,δds,color=δds,markersize = 0.1)
# fig