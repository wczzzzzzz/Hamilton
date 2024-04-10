
using Revise, ApproxOperator, Printf, SparseArrays

include("import_Scordelis_Lo_roof.jl")

ndiv= 11
elements,nodes = import_roof_Tri3("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

ρA = 1
EA = 1
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)


ops = [
    Operator{:∫qmpdΩ}(:ρA=>ρA),
    Operator{:∫qkpdΩ}(:EA=>EA),
]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)

ops[1](elements["Ω"],m)
ops[2](elements["Ω"],k)

Δt = 5
d = zeros(nₚ,nₜ)
d̈ₙ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

for n in 1:nₚ
    global d̈ₙ .+= m/k *d[:,n] 
    global ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ
    global d[:,n+1] .= d[:,n] + Δt*ḋₙ

    # for i in (1:nₚ)
    # global d₁₁ₙ .+= m/k *d[:, i] 
    # global d₁ₙ₊₁ .+= d₁ₙ + Δt*d₁₁ₙ
    # global dₙ₊₁ .+= dₙ + Δt*d₁ₙ
end
