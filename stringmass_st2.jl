using CairoMakie, LinearAlgebra

𝑘 = 100
𝑚 = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.005:8.0
𝜔 = (𝑘/𝑚)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

t = 0.0:0.005:8.0
nₚ = length(t)
nₑ = nₚ-1

# FEM
k_ = zeros(nₚ,nₚ)

for i in 1:nₑ
    t₁ = t[i]
    t₂ = t[i+1]
    𝐿 = t₂ - t₁
    k_[i,i] += 𝑚/𝐿 - 𝑘/3*𝐿
    k_[i,i+1] += -𝑚/𝐿 - 𝑘/6*𝐿
    k_[i+1,i] += -𝑚/𝐿 - 𝑘/6*𝐿
    k_[i+1,i+1] += 𝑚/𝐿 - 𝑘/3*𝐿
end

# F1 formulation: δq₀ = 0, δq̇₀ = 0
k = zeros(nₚ,nₚ)
k .= k_
f = zeros(nₚ)
𝐿 = t[nₚ]-t[nₚ-1]
k[nₚ,nₚ-1] -= -𝑚/𝐿
k[nₚ,nₚ] -= 𝑚/𝐿

𝐿 = t[2]-t[1]
f .-= k[:,1]*q₀
f .-= k[:,2]*(q₀+𝐿*q̇₀)
k[1,:] .= 0.0
k[2,:] .= 0.0
k[:,1] .= 0.0
k[:,2] .= 0.0
k[1,1] = 1.0
k[2,2] = 1.0
f[1] = q₀
f[2] = q₀+𝐿*q̇₀

d = k\f

lines!(t, d, color = :blue)

# F2 formulation: δq₀ = 0, δq₁ = 0
k = zeros(nₚ,nₚ)
k .= k_
f = zeros(nₚ)
𝐿 = t[2]-t[1]
f .-= k[:,1]*q₀
f .-= k[:,2]*(q₀+𝐿*q̇₀)
k[1,:] .= 0.0
k[nₚ,:] .= 0.0
k[:,1] .= 0.0
k[:,2] .= 0.0
k[1,1] = 1.0
k[nₚ,2] = 1.0
f[1] = q₀
f[nₚ] = q₀+𝐿*q̇₀

d = k\f

# lines!(t, d, color = :orange)

fig

# val = eigvals(k)