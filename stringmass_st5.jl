using  ApproxOperator, JuMP, Ipopt, CairoMakie

model = Model(Ipopt.Optimizer)

include("import_Scordelis_Lo_roof.jl")

ndiv= 11
elements,nodes = import_roof_Tri3("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

α = 1e9
ρA = 1
EA = 1
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->𝑇(y))
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->𝑇(y))


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
ops[2](elements["Γ₁"],f)
ops[3](elements["Γ₄"],f)
ops[4](elements["Γ₁"],kᵅ,fᵅ)
ops[4](elements["Γ₂"],kᵅ,fᵅ)
ops[4](elements["Γ₃"],kᵝ,fᵝ)

d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
# δd = d[nₚ+1:end]
# d = d[1:nₚ]

push!(nodes,:d=>d)
fig = Figure()
Axis(fig[1, 1])
xs = [node.x for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
ys = [node.d for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
lines!(xs,ys, color = :blue)



fig