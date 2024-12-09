using  ApproxOperator, CairoMakie

# model = Model(Ipopt.Optimizer)

include("import_hmd.jl")

ndiv= 160
elements,nodes = import_hmd_bar("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

kᶜ = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.05:8.0
𝜔 = (kᶜ/m)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)


       Operator{:∫q̇mpqkpdx}(:m=>m,:kᶜ=>kᶜ),
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

# ops[2](elements["Γᵍ"],k,f)

ops[1](elements["Ω"],k)

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

α = 1e9
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
kα[1,1] += α
fα[1] += α*q₀
kβ = zeros(nₚ,nₚ)
kβ[nₚ,nₚ] += α

d = [k+kα -k;-k kβ]\[fα;f]
δd = d[nₚ+1:end]
d = d[1:nₚ]


lines!(nodes.x, d, color = :blue)

# e = d - 𝑥
# lines!(𝑡, e, color = :red)

fig

# save("./fig/一维/hmd_1d.png",fig)
