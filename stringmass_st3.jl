using CairoMakie

𝑘 = 1e3
𝑚 = 1.0
q̇₀ = 5.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.01:1.0
𝜔 = (𝑘/𝑚)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

dt = 0.1
t = collect(0.0:dt:1.0)
nₚ = length(t)
nₑ = nₚ-1
# for i in 2:nₚ-1
#     t[i] += dt*(rand()-0.5)
# end

# FEM weak LM
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

α = 1e8
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*q₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵝ[nₚ,nₚ] += α

d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
d = d[1:nₚ]

lines!(t, d, color = :blue)

# FEM weak test
# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)

# for i in 1:nₑ
#     t₁ = t[i]
#     t₂ = t[i+1]
#     𝐿 = t₂ - t₁
#     k[i,i] += 𝑚/𝐿 - 𝑘/3*𝐿
#     k[i,i+1] += -𝑚/𝐿 - 𝑘/6*𝐿
#     k[i+1,i] += -𝑚/𝐿 - 𝑘/6*𝐿
#     k[i+1,i+1] += 𝑚/𝐿 - 𝑘/3*𝐿
# end

# kᵥ = k[1:nₚ-1,1:nₚ-1]
# fᵥ = f[1:nₚ-1]

# 𝑃₀ = 𝑚*q̇₀
# f[1] -= 𝑃₀
# α = 1e8
# k[1,1] += α
# f[1] += α*q₀

# d = k\f
# d[1:nₚ-1] .-=  kᵥ\fᵥ

# lines!(t, d, color = :red)

fig