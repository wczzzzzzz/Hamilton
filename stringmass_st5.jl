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
prescribe!(elements["Γᵍ"],:𝑃=>(x,y,z)->100.0)
# prescribe!(elements["Ω"],:b=>(x,y,z)->0.0)

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.01:1.0
𝜔 = (kᶜ/m)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops = [
       
       Operator{:∫∫q̇mpqkpdx}(:ρA=>ρA,:EA=>EA),
       Operator{:∫𝑃δudx}(),
]



ops[1](elements["Ω"],k)
ops[2](elements["Γᵍ"],f)

α = 1e9
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)
kα[1,1] += α
fα[1] += α*q₀
kβ = zeros(nₚ,nₚ)
kβ[nₚ,nₚ] += α

d = [k+kα k;k kβ]\[f+fα;f]
δd = d[nₚ+1:end]
d = d[1:nₚ]


lines!(nodes.x, d, color = :blue)
# lines!(t, d, color = :blue)


fig