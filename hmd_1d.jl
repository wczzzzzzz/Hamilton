using  ApproxOperator, CairoMakie
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫q̇mṗqkpdxdt

# model = Model(Ipopt.Optimizer)

include("import_hmd.jl")

ndiv= 160
elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

kᶜ = 1e2
m = 1.0
q̇₀ = 5.0
q₀ = 1.0
prescribe!(elements["Ω"],:m=>(x,y,z)->m)
prescribe!(elements["Ω"],:kᶜ=>(x,y,z)->kᶜ)

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.005:8.0
𝜔 = (kᶜ/m)^0.5
# 𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
# lines!(𝑡, 𝑥, color = :black)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

𝑎 = ∫∫q̇mṗqkpdxdt=>elements["Ω"]

𝑎(k)

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

α = 1e12
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*q₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
# kᵝ[1,1] += α
kᵝ[nₚ,nₚ] += α

d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = d[1:nₚ]


lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2]], color = :blue)
# lines!(nodes.x, d, color = :blue)

# e = d - 𝑢.(t)
# e = d - 𝑥

# lines!(𝑡, e[[1,3:end...,2]], color = :red)


fig

# save("./fig/一维/hmd_1d.png",fig)
