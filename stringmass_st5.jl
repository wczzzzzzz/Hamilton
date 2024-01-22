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

# kᶜ = 100.0
# m = 1.0
# q̇₀ = 5.0
# q₀ = 1.0
ρ = 1.0
A = 100.0
E = 3e6
EA = E*A
# prescribe!(elements["Γᵗ"],:P=>(x,y,z)->1.0)
# prescribe!(elements["Ω"],:b=>(x,y,z)->0.0)

# fig = Figure()
# Axis(fig[1, 1])
# 𝑡 = 0.0:0.01:1.0
# 𝜔 = (kᶜ/m)^0.5
# 𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
# lines!(𝑡, 𝑥, color = :black)


ops = [
       Operator{:∫∫q̇mpqkpdx}(:ρ=>ρ,:A=>A,:EA=>EA),
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k,f)
# ops[2](elements["Γᵍ"],k,f)

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


lines!(nodes.x[[1,3:end...,2]], d, color = :blue)


fig