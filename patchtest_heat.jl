using  ApproxOperator, JuMP, Ipopt, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy
import ApproxOperator.Test: cc𝝭, cc∇𝝭
using GLMakie
using SparseArrays

include("import_hmd.jl")

ndiv= 32
# elements,nodes,nodes_t = import_hermite("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
elements,nodes = import_hmd_Tri3("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"

nₚ = length(nodes)
# nₑ = length(elements["Ωᵗ"])
# nₗ = length(nodes_t) - nₚ - nₑ

# set𝝭!(elements["Ωᵗ"])
# set∇𝝭!(elements["Ωᵗ"])
# set𝝭!(elements["Γ₁ᵗ"])
# set𝝭!(elements["Γ₂ᵗ"])
# set𝝭!(elements["Γ₃ᵗ"])
# set𝝭!(elements["Γ₄ᵗ"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

α = 1e19
n = 3
u(x,y) = (x+y)^n
v(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

# prescribe!(elements["Ωᵗ"],:k=>(x,y,z)->1.0,index=:𝑔)
# prescribe!(elements["Ωᵗ"],:b=>b)
# prescribe!(elements["Ωᵗ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵗ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵗ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# prescribe!(elements["Γ₄ᵗ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Γ₂ᵗ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Γ₁ᵗ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Γ₁ᵗ"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₂ᵗ"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₃ᵗ"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₄ᵗ"],:α=>(x,y,z)->α)

prescribe!(elements["Ω"],:k=>(x,y,z)->1.0,index=:𝑔)
prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Γ₄"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ₃"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ₂"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ₁"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)

# k = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# f = zeros(nₚ+nₗ+nₑ)
# fᵦ = zeros(nₚ+nₗ+nₑ)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fᵦ = zeros(nₚ)

# 𝑎 = ∫∫∇v∇udxdy=>elements["Ωᵗ"]
# 𝑓ᵦ = ∫vbdΩ=>elements["Ωᵗ"]
# 𝑓 = ∫vgdΓ=>elements["Γ₂ᵗ"]∪elements["Γ₃ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₁ᵗ"]

𝑎 = ∫∫∇v∇udxdy=>elements["Ω"]
𝑓ᵦ  = ∫vbdΩ=>elements["Ω"]
𝑓 = ∫vgdΓ=>elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₁"]

𝑎(k)
𝑓(k,f)
𝑓ᵦ(fᵦ)

d=k\(f+fᵦ)

push!(nodes,:d=>d)
# push!(nodes_t,:d=>d)

𝐿₂ = log10.(L₂(elements["Ω"]))
# 𝐿₂ = log10.(L₂(elements["Ωᵗ"]))


index = [4,8,16,32]
# # index = [0.4,0.3,0.2,0.1]
# # index = [0,1,2,3]
XLSX.openxlsx("./excel/patchtest.xlsx", mode="rw") do xf
    Sheet = xf[1]
    ind = findfirst(n->n==ndiv,index)+1
    Sheet["A"*string(ind)] = log10(4/ndiv)
    # Sheet["A"*string(ind)] = log10(nₚ)
    Sheet["B"*string(ind)] = 𝐿₂
end