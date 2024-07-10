using CairoMakie, LinearAlgebra

𝑘 = 100.0
𝑚 = 1.0
q̇₀ = 5.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.001:1.0
𝜔 = (𝑘/𝑚)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
u(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
lines!(𝑡, 𝑥, color = :black)

t = 0.0:0.001:1.0
nₚ = length(t)
nₑ = nₚ-1

# FEM weak LM
kᵤᵤ = zeros(nₚ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
kᵤₚ = zeros(nₚ,nₚ)
fᵤ = zeros(nₚ)
fₚ = zeros(nₚ)

for i in 1:nₑ
    t₁ = t[i]
    t₂ = t[i+1]
    𝐿 = t₂ - t₁

    kᵤᵤ[i,i] += 𝐿*𝑘/3
    kᵤᵤ[i,i+1] += 𝐿*𝑘/6
    kᵤᵤ[i+1,i] += 𝐿*𝑘/6
    kᵤᵤ[i+1,i+1] += 𝐿*𝑘/3
    
    kₚₚ[i,i] += 𝐿/𝑚/3
    kₚₚ[i,i+1] += 𝐿/𝑚/6
    kₚₚ[i+1,i] += 𝐿/𝑚/6
    kₚₚ[i+1,i+1] += 𝐿/𝑚/3

    kᵤₚ[i,i] -= -0.5
    kᵤₚ[i,i+1] -= -0.5
    kᵤₚ[i+1,i] -= 0.5
    kᵤₚ[i+1,i+1] -= 0.5
end

# α = 1e9
# kᵤᵤ[1,1] += α
# fᵤ[1] += α*q₀
# kₚₚ[1,1] += α
# fₚ[1] += α*𝑚*q̇₀
# kᵤᵤ[1,1] += α
# fᵤ[1] += α*𝑚*q̇₀

kᵤₚ[1,1] = 1.5
kᵤₚ[nₚ,nₚ] = 0.5
fₚ[1] = -1
fₚ[nₚ] = cos(10) + 0.5*sin(10)
# fₚ[nₚ] = cos(10/180*π) + 0.5*sin(10/180*π)


k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [fᵤ;fₚ]

d = k\f

lines!(t, d[1:nₚ], color = :blue)

fig

# val = eigvals(kᵤᵤ)
# val = eigvals(kₚₚ)
val = eigvals(kᵤₚ)