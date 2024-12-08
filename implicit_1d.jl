using   GLMakie, LinearAlgebra

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "T", ylabel = "x",title = "implicit Euler vs Exact Solution")
𝑡 = 0.0:0.05:8.0
𝜔 = (k/m)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(ax, 𝑡, 𝑥, color = :black)

T = 8.0
Δt = 0.05
t = 0.0:Δt:T
nₜ = Int(T/Δt)

q = zeros(nₜ+1)
q̇ = zeros(nₜ+1)
q̈ = zeros(nₜ+1)

q[1] = q₀
q̇[1] = q̇₀

for n in 1:nₜ
    q̈[n+1] = -k/m * q[n+1]  
    q̇[n+1] = q̇[n] + Δt * q̈[n+1]  
    q[n+1] = (m*q[n])/(m + k*Δt^2) + (m*Δt*q̇[n])/(m + k*Δt^2)
end

lines!(ax, t, q, color = :blue)

# e = q - 𝑥
# lines!(t, e, color = :red)
xlims!(ax, 0, 8)

fig

# save("./fig/一维/semi_8.png",fig)
# save("./fig/一维/e.png",fig)

