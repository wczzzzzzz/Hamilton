using CairoMakie

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
q̇ = zeros(nₜ+1)
c = zeros(nₜ+1)
q[1] = q₀
c[1] = m * q̇₀  

# for n in 1:nₜ
#     c[n+1] = c[n] - k * q[n] * Δt  
#     q[n+1] = q[n] + c[n+1] / m * Δt  
# end

for n in 1:2*nₜ
    q[n+1] = q[n] + c[1] / m * Δt
    c[n+2] = c[n] - k * q[n+1] * Δt  
    q[n+2] = q[n+1] + c[n+2] / m * Δt  
en

# e = q - 𝑥
# lines!(t, e, color = :red)

lines!(ax, t, q, color = :blue)
xlims!(ax, 0, 8)

fig