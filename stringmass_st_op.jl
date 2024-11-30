using CairoMakie, LinearAlgebra, JuMP

# using HiGHS
# model = Model(HiGHS.Optimizer)

using Ipopt
model = Model(Ipopt.Optimizer)

𝑘 = 100.0
𝑚 = 1.0
q̇₀ = 5.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.01:1.0
𝜔 = (𝑘/𝑚)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

t = 0.0:0.01:1.0
nₚ = length(t)
nₑ = nₚ-1

# FEM optimization
k = zeros(nₚ,nₚ)
f = zeros(nₚ)

for i in 1:nₑ
    t₁ = t[i]
    t₂ = t[i+1]
    𝐿 = t₂ - t₁
    k[i,i] += 𝑚/𝐿 - 𝑘/3*𝐿
    k[i,i+1] += -𝑚/𝐿 - 𝑘/6*𝐿
    k[i+1,i] += -𝑚/𝐿 - 𝑘/6*𝐿
    k[i+1,i+1] += 𝑚/𝐿 - 𝑘/3*𝐿
end

𝑃₀ = 𝑚*q̇₀
f[1] -= 𝑃₀
kc = k[1:nₚ-1,1:nₚ]
fc = f[1:nₚ-1]
𝐿 = t[2]-t[1]

# Optimization
@variable(model,d[1:nₚ])
@constraint(model, kc*d .== fc)
@constraint(model, d[1] == q₀)
# @constraint(model, (d[2]-d[1])/𝐿 == q̇₀)
@objective(model, Min, d'*(0.5.*k*d-f))
optimize!(model)
d = JuMP.value.(d)
lines!(t, d, color = :blue)

# FEM
# 𝑃₀ = 𝑚*q̇₀
# f[1] -= 𝑃₀
# f .-= k[:,1].*q₀
# k[:,1] .= 0.0
# k[nₚ,:] .= 0.0
# k[nₚ,1] = 1.0
# f[nₚ] = q₀

# d = k\f

# lines!(t, d, color = :orange)

fig