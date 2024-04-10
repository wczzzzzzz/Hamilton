
using Revise, ApproxOperator, Printf, SparseArrays

include("import_hmd_test.jl")

ndiv= 4
elements,nodes = import_hmd_bar("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₜ = length(nodes)

set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

ρA = 1
EA = 1
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:g=>(x,y,z)->0.0)

ops = [
    Operator{:∫qmpdΩ}(:ρA=>ρA),
    Operator{:∫qkpdΩ}(:EA=>EA),
]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)

ops[1](elements["Ω"],m)
ops[2](elements["Ω"],k)

Δt = 5
d = zeros(nₚ,nₜ+1)
d̈ₙ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

for n in 1:nₚ
    global d̈ₙ .+= k/m *d[:,n] 
    global ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ
    global d[:,n+1] .= d[:,n] + Δt*ḋₙ

    # for i in (1:nₚ)
    # global d₁₁ₙ .+= m/k *d[:, i] 
    # global d₁ₙ₊₁ .+= d₁ₙ + Δt*d₁₁ₙ
    # global dₙ₊₁ .+= dₙ + Δt*d₁ₙ
end
