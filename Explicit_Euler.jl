
using Revise, ApproxOperator, Printf, SparseArrays, LinearAlgebra, CairoMakie

include("import_hmd.jl")

ndiv= 80
elements,nodes = import_hmd_bar("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

α = 1e9
ρA = 1
EA = 1
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)

fig = Figure()
Axis(fig[1, 1])

ops = [
    Operator{:∫qmpdΩ}(:ρA=>ρA),
    Operator{:∫qkpdΩ}(:EA=>EA),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>α),
]

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
fᵗ = zeros(nₚ)
fᵍ = zeros(nₚ)

ops[1](elements["Ω"],m)
ops[2](elements["Ω"],k)
ops[4](elements["Γᵍ"],m,fᵍ)

T = 4
Δt = 0.1
nₜ = Int(T/Δt)
d = zeros(nₚ,nₜ+1)
d̈ₙ = zeros(nₚ)
ḋₙ = zeros(nₚ)

𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)


for n in 1:nₜ
    fill!(fᵗ,0.0)
    t = n*Δt
    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->𝑇(t))
    ops[3](elements["Γᵗ"],fᵗ)

     d̈ₙ .= m\(fᵗ+fᵍ - k*d[:,n])
     d[:,n+1] .= d[:,n] + Δt*ḋₙ
     ḋₙ .+= Δt*d̈ₙ
end

lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2],21], color = :blue)

fig