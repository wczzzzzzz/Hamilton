using CairoMakie, LinearAlgebra

𝑘 = 100
𝑚 = 1.0
q̇₀ = 5.0
q₀ = 1.0

fig = Figure()
# Axis(fig[1, 1])
ax = Axis(fig[1, 1], xlabel = "T", ylabel = "x",title = "Hamilton vs Exact Solution")
𝑡 = 0.0:0.05:8.0
𝜔 = (𝑘/𝑚)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
𝑝 = 𝑚*𝜔.*(q̇₀/𝜔.*cos.(𝜔.*𝑡)-q₀.*sin.(𝜔.*𝑡))
u(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
p(t) = 𝑚*𝜔*(q̇₀/𝜔*cos(𝜔*t) - q₀*sin(𝜔*t))
lines!(ax, 𝑡, 𝑥, color = :black)

# lines!(𝑡, 𝑝, color = :black)

t = 0.0:0.05:8.0
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

kᵤₚ[1,1] += -1.0
fₚ[1] += -u(0.)
kᵤₚ[nₚ,nₚ] += 1.0
fₚ[nₚ] += u(1.)
# fₚ[nₚ] += u(1.)
# fᵤ[1] +=  p(0.)
# fᵤ[nₚ] += -p(1.)
# kₚₚ[nₚ,nₚ] += 1e9

k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [fᵤ;fₚ]
# k[nₚ,2*nₚ] += 1.
# k[nₚ+1,1] += 1.
# k[2*nₚ,nₚ] += 1.
# f[nₚ+1] += u(0.)
# f[2*nₚ] += u(1.)

# kₗ = zeros(nₚ)
# kₗ[1,1] = 1.
# fₗ = [u(0.)]
# k = [kᵤᵤ kᵤₚ kₗ;kᵤₚ' kₚₚ zeros(nₚ);kₗ' zeros(nₚ)' [0.]]
# f = [fᵤ;fₚ;fₗ]

d = k\f

# val = eigvals(kᵤᵤ)
# val = eigvals(kₚₚ)
# val = eigvals(kᵤₚ)
val = eigvals(k)
# val = eigvals(kᵤₚ*inv(kₚₚ)*kᵤₚ')
e = d[1:nₚ] - 𝑥
# lines!(ax, t, e, color = :red)
lines!(ax, t, d[1:nₚ], color = :blue)

# lines!(t, d[nₚ+1:end-1], color = :blue)

# invisible_line = lines!(ax, [0, 0], [0, 0], color = :white, label="Δt=0.05", visible=false)
# blue_line = lines!(ax, t[1:2], q[1:2], color = :blue, label="Hamilton")
# black_line = lines!(ax, t[1:2], 𝑢.(t[1:2]), color = :black, label="Exact Solution")
# red_line = lines!(ax, t[1:2], e[1:2], color = :red, label="error")
# leg = Legend(fig[1, 2], [red_line, invisible_line], ["error", "Δt=0.05"], position=(0.95, 0.95))
# leg = Legend(fig[1, 2], [blue_line, black_line, invisible_line], ["Hamilton", "Exact Solution", "Δt=0.05"], position=(0.95, 0.95))

fig

# save("./fig/一维/hmd_1d.png",fig)
# save("./fig/一维/hmd_1d_e.png",fig)