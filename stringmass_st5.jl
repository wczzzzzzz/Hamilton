using  ApproxOperator, JuMP, Ipopt, CairoMakie

model = Model(Ipopt.Optimizer)

include("import_Scordelis_Lo_roof.jl")

ndiv= 11
elements,nodes = import_roof_Tri3("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

kᶜ = 100.0
m = 1.0
q̇₀ = 5.0
q₀ = 1.0
ρ = 1.0
A = 100.0
E = 3e6
ρA = ρ*A
EA = E*A
𝑇(t) = t > 1.0 ? 0 : - sin(π*t)
prescribe!(elements["Γᵍ"],:𝑃=>(x,y,z)->100.0)
prescribe!(elements["Γ₄"],:𝑃=>(x,y,z)->𝑇(y))
# prescribe!(elements["Ω"],:b=>(x,y,z)->0.0)

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.01:1.0
𝜔 = (kᶜ/m)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

ops = [
       
       Operator{:∫∫q̇mpqkpdx}(:ρA=>ρA,:EA=>EA),
       Operator{:∫𝑃δudx}(),
       Operator{:∫vtdΓ}(),
       Operator{:∫vgdΓ}(:α=>α),
]



ops[1](elements["Ω"],k)
ops[2](elements["Γᵍ"],f)
ops[3](elements["Γ₄"],f)
ops[4](elements["Γ₁"],kᵅ,fᵅ)
ops[4](elements["Γ₂"],kᵅ,fᵅ)
ops[4](elements["Γ₃"],kᵝ,fᵝ)

d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]

δd = d[nₚ+1:end]
d = d[1:nₚ]


lines!(nodes.x, d, color = :blue)
# lines!(t, d, color = :blue)


fig