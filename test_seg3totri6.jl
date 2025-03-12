
using ApproxOperator
import ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫vtdΓ, ∫∇𝑛vgds
import Gmsh: gmsh

filename = "msh/tri6_x=20/20.msh"
gmsh.initialize()
gmsh.open(filename)
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"])
elements["Γ₁"] = getElements(nodes, entities["Γ¹"], normal = true)
elements["Γ₂"] = getElements(nodes, entities["Γ²"], normal = true)
elements["Γ₃"] = getElements(nodes, entities["Γ³"], normal = true)
elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], normal = true)
elements["Ω∩Γ₁"] = ApproxOperator.Seg3toTri6(elements["Γ₁"],elements["Ω"])
elements["Ω∩Γ₂"] = ApproxOperator.Seg3toTri6(elements["Γ₂"],elements["Ω"])
elements["Ω∩Γ₃"] = ApproxOperator.Seg3toTri6(elements["Γ₃"],elements["Ω"])
elements["Ω∩Γ₄"] = ApproxOperator.Seg3toTri6(elements["Γ₄"],elements["Ω"])

nₚ = length(nodes)
push!(elements["Ω"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
push!(elements["Γ₁"], :𝝭)
push!(elements["Γ₂"], :𝝭)
push!(elements["Γ₃"], :𝝭)
push!(elements["Γ₄"], :𝝭)
push!(elements["Ω∩Γ₁"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
push!(elements["Ω∩Γ₂"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
push!(elements["Ω∩Γ₃"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
push!(elements["Ω∩Γ₄"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ω∩Γ₁"])
set∇𝝭!(elements["Ω∩Γ₂"])
set∇𝝭!(elements["Ω∩Γ₃"])
set∇𝝭!(elements["Ω∩Γ₄"])

u(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0

prescribe!(elements["Ω"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ₁"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Γ₂"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Γ₃"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Γ₄"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Ω∩Γ₁"],:α=>(x,y,z)->0.0)
prescribe!(elements["Ω∩Γ₂"],:α=>(x,y,z)->0.0)
prescribe!(elements["Ω∩Γ₃"],:α=>(x,y,z)->0.0)
prescribe!(elements["Ω∩Γ₄"],:α=>(x,y,z)->0.0)
prescribe!(elements["Ω∩Γ₁"],:k=>(x,y,z)->1.0)
prescribe!(elements["Ω∩Γ₂"],:k=>(x,y,z)->1.0)
prescribe!(elements["Ω∩Γ₃"],:k=>(x,y,z)->1.0)
prescribe!(elements["Ω∩Γ₄"],:k=>(x,y,z)->1.0)
prescribe!(elements["Ω∩Γ₁"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ω∩Γ₂"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ω∩Γ₃"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ω∩Γ₄"],:g=>(x,y,z)->u(x,y))

𝑎 = ∫∫∇v∇udxdy=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
# 𝑓ₙ = ∫∇𝑛vgds=>elements["Ω∩Γ₁"][20:20]∪elements["Ω∩Γ₂"][1:1]
# 𝑓ₙ = ∫∇𝑛vgds=>elements["Ω∩Γ₁"][20:20]
# 𝑓ₙ = ∫∇𝑛vgds=>elements["Ω∩Γ₄"]
# 𝑓ₙ = ∫∇𝑛vgds=>elements["Ω∩Γ₂"]
𝑓ₙ = ∫∇𝑛vgds=>elements["Ω∩Γ₁"]∪elements["Ω∩Γ₂"]∪elements["Ω∩Γ₃"]∪elements["Ω∩Γ₄"]

k = zeros(nₚ,nₚ)
kₙ = zeros(nₚ,nₚ)
f = zeros(nₚ)
fₙ = zeros(nₚ)

𝑎(k)
𝑓(f)
d = [u(node.x,node.y) for node in nodes]

err1 = k*d - f
𝑓ₙ(kₙ,fₙ)
err2 = (k+kₙ)*d - fₙ