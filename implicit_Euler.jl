
using Revise, ApproxOperator, Printf, SparseArrays, LinearAlgebra, CairoMakie, XLSX

include("import_hmd.jl")

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

α = (EA/ρA)^0.5
function 𝑢(x,t)
    if x < α*(t-1)
        return 2*α/π
    elseif α*t < x
        return 0
    else
        α/π*(1-cos(π*(t-x/α)))
    end
end



for n in 1:nₜ
    fill!(fᵗ,0.0)
    t = (n+1)*Δt
    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->-𝑇(t))
    ops[3](elements["Γᵗ"],fᵗ)

     d̈ₙ₊₁ .= m\(fᵗ+fᵍ - k*d[:,n+1])
     ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ₊₁
     d[:,n+1] .= (m + k*Δt^2)\m*d[:,n] + (m + k*Δt^2)\(m*Δt)*ḋₙ + (m + k*Δt^2)/(Δt^2)*(fᵗ+fᵍ)
end

ys = 0.0:4.0/(41-1):4.0

# for (i, node) in enumerate(nodes)
#     for (j, t) in enumerate(ys)
#         x = node.x
#         z = d[i,j]
#         Δ = d[i,j] - 𝑢(x,t)
#         XLSX.openxlsx("./excel/implicit_Euler.xlsx", mode="rw") do xf
#             Sheet = xf[1]
#             ind = findfirst(n->n==ndiv,20)+(i-1)*41+j
#             Sheet["A"*string(ind)] = x
#             Sheet["B"*string(ind)] = t
#             Sheet["C"*string(ind)] = z
#             Sheet["D"*string(ind)] = Δ
#         end
#     end
# end

lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2],21], color = :blue)

fig