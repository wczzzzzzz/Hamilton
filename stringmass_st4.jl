using  ApproxOperator, JuMP, Ipopt

model = Model(Ipopt.Optimizer)

include("import_Scordelis_Lo_roof.jl")

ndiv= 3
elements,nodes = import_roof_gauss("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

kᶜ = 100.0
m = 1.0
q̇₀ = 5.0
q₀ = 1.0

ops = [
       Operator{:∫q̇mpqkpdx}(:m=>m,:kᶜ=>kᶜ),
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

# ops[2](elements["Γᵍ"],k,f)

ops[1](elements["Ω"],k)

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

d = k\f


lines!(t, d, color = :blue)


