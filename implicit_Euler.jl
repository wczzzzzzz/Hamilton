
using Revise, ApproxOperator, Printf, SparseArrays, LinearAlgebra, CairoMakie

include("import_hmd_test.jl")

ndiv= 20
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
d̈ₙ₊₁ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)


for n in 1:nₜ
    fill!(fᵗ,0.0)
    t = (n+1)*Δt
    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->-𝑇(t))
    ops[3](elements["Γᵗ"],fᵗ)

     d̈ₙ₊₁ .= m\(fᵗ+fᵍ - k*d[:,n+1])
     ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ₊₁
     d[:,n+1] .= (m + k*Δt^2)\m*d[:,n] + (m + k*Δt^2)\(m*Δt)*ḋₙ + (m + k*Δt^2)/(Δt^2)*(fᵗ+fᵍ)

    #  XLSX.openxlsx("./excel/implicit_Euler.xlsx", mode="rw") do xf
    #     Sheet = xf[1]
    #     ind = findfirst(n->n==ndiv,20)+1
    #     Sheet["B"*string(ind)] = d
    # end
end

# for i in 1:21
#     x = nodes.x[i]
#     y = nodes.y[i]
#          XLSX.openxlsx("./excel/implicit_Euler.xlsx", mode="rw") do xf
#         Sheet = xf[2]
#         ind = findfirst(n->n==ndiv,20)+i
#             Sheet["C"*string(ind)] = x
#             Sheet["D"*string(ind)] = y
        
#     end
# end

lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2],21], color = :blue)

fig