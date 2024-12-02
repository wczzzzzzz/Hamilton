using GLMakie

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "t", ylabel = "x", title = "Leapfrog Method vs Exact Solution")
𝑡 = 0.0:0.05:8.0
𝜔 = (k/m)^0.5
𝑥 = q₀ .* cos.(𝜔 .* 𝑡) + q̇₀ / 𝜔 .* sin.(𝜔 .* 𝑡)
lines!(ax, 𝑡, 𝑥, color = :black)

T = 8.0
Δt = 0.05
t = 0.0:Δt:T
nₜ = Int(T/Δt)

q = zeros(nₜ+1)
q[1] = q₀
q̇ = zeros(nₜ+1)
q̇[1] = q̇₀

points = Observable(Point2f[(0, q[1])])

scatter!(ax, points, color = :blue)

scatter!(ax, points, color = :blue, markersize=10)
record(fig, "./fig/一维/dian.mp4", 1:nₜ) do frame
    q̈ = -k/m * q[frame]
    q̇[frame+1] = q̇[frame] + Δt * q̈
    q[frame+1] = q[frame] + Δt * q̇[frame+1]
    push!(points[], Point2f(𝑡[frame], q[frame+1]))
    scatter!(ax, points, color = :blue)
end


# record(fig, "./fig/一维/xian.mp4", 1:nₜ) do frame
#     q̈ = -k/m * q[frame]
#     q̇[frame+1] = q̇[frame] + Δt * q̈
#     q[frame+1] = q[frame] + Δt * q̇[frame+1]
#     push!(points[], Point2f(𝑡[frame], q[frame+1]))
#     lines!(ax, t, q, color = :blue)  
# end

# qx = Observable([t])
# qy = Observable([q])

# lines!(ax, qx, qy, color = :blue)

# record(fig, "spring_mass_simulation.mp4", 1:nₜ) do frame
#     q̈ = -k/m * q[frame]
#     q̇[frame+1] = q̇[frame] + Δt * q̈
#     q[frame+1] = q[frame] + Δt * q̇[frame+1]
#     push!(points[], Point2f(qx[], t[frame]))
#     push!(points[], Point2f(qy[], q[frame+1]))
# end

xlims!(ax, 0, 8)
fig