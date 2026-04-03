using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, Hₑ
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ

using TimerOutputs, WriteVTK
# using LinearSolve
import Gmsh: gmsh

EA = 1.0
ρA = 1.0
b = 1.0
c = (EA/ρA)^0.5/b
𝑢(x,t) = sin(x-c*t)
𝑢̇(x,t) = -c*cos(x-c*t)
𝑢ₓ(x,t) = cos(x-c*t)

δ𝑢(x,t) = 0.0
δ𝑢̇(x,t) = 0.0
δ𝑢ₓ(x,t) = 0.0

# δ𝑢(x,t) = sin(x-c*t)
# δ𝑢̇(x,t) = -c*cos(x-c*t)
# δ𝑢ₓ(x,t) = cos(x-c*t)

# δ𝑢(x,t) = sin(x)*cos(c*t)
# δ𝑢̇(x,t) = -c*sin(x)*sin(c*t)
# δ𝑢ₓ(x,t) = cos(x)*cos(c*t)

const to = TimerOutput()
gmsh.initialize()

# nˣ = 10
# nʸ = 13
# for nʸ in 4:16
filename = "square_4_refined_r3"
# filename = "square_tri3_irregular_1_r3"
# filename = "square_tri3_8"
@timeit to "open msh file" gmsh.open("msh/$filename.msh")
# @timeit to "open msh file" gmsh.open("msh/square_tri3_$n.msh")
# @timeit to "open msh file" gmsh.open("msh/square_tri3_$nˣ"*"_$nʸ.msh")
# @timeit to "open msh file" gmsh.open("msh/square_tri6_$nˣ"*"_$nʸ.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

integrationOrder = 2
@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" element = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(element, :EA=>EA, :ρA=>ρA)
    @timeit to "calculate shape functions" set∇𝝭!(element)
    𝑎 = ∫∫∇q∇pdxdt=>element
    @timeit to "assemble" 𝑎(k)
end

α = 1e8
@timeit to "calculate boundary for u" begin
    @timeit to "get elements" element_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" element_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" element_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" element_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(element_1,:g=>(x,y,z)->𝑢(x,y), :α=>α)
    prescribe!(element_2,:g=>(x,y,z)->𝑢(x,y), :α=>α)
    prescribe!(element_3,:g=>(x,y,z)->𝑢(x,y), :α=>α)
    prescribe!(element_4,:g=>(x,y,z)->𝑢(x,y), :α=>α)
    @timeit to "calculate shape functions" set𝝭!(element_1)
    @timeit to "calculate shape functions" set𝝭!(element_2)
    @timeit to "calculate shape functions" set𝝭!(element_3)
    @timeit to "calculate shape functions" set𝝭!(element_4)
    𝑎ᵅ = ∫vgdΓ=>element_1∪element_2∪element_4
    # 𝑎ᵅ = ∫vgdΓ=>element_1∪element_2∪element_3∪element_4
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
end
@timeit to "calculate δu" begin
    @timeit to "get elements" element_1 = getElements(nodes, entities["Γ¹"],integrationOrder,normal=true)
    @timeit to "get elements" element_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" element_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" element_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(element_1,:t=>(x,y,z,n₁,n₂)->ρA*𝑢̇(x,y)*n₂)
    prescribe!(element_1,:g=>(x,y,z)->δ𝑢(x,y), :α=>α)
    prescribe!(element_2,:g=>(x,y,z)->δ𝑢(x,y), :α=>α)
    prescribe!(element_3,:g=>(x,y,z)->δ𝑢(x,y), :α=>α)
    prescribe!(element_4,:g=>(x,y,z)->δ𝑢(x,y), :α=>α)
    @timeit to "calculate shape functions" set𝝭!(element_1)
    @timeit to "calculate shape functions" set𝝭!(element_2)
    @timeit to "calculate shape functions" set𝝭!(element_3)
    @timeit to "calculate shape functions" set𝝭!(element_4)
    𝑓 = ∫vtdΓ=>element_1
    𝑎ᵝ = ∫vgdΓ=>element_2∪element_3∪element_4
    # 𝑎ᵝ = ∫vgdΓ=>element_1∪element_2∪element_3∪element_4
    @timeit to "assemble" 𝑓(f)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

# dt = (k+kᵅ)\fᵅ
# d = dt
# δd = zeros(nₚ)

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
δd = dt[nₚ+1:end]

# prob = LinearProblem(k+kᵅ,fᵅ)
# prob = LinearProblem([k+kᵅ -k;-k kᵝ],[fᵅ;-f+fᵝ])
# linsolve = init(prob)
# sol = solve!(linsolve)

# d = sol.u[1:nₚ]
# δd = sol.u[nₚ+1:end]
# δd = zeros(nₚ)

push!(nodes, :d=>d)
push!(nodes,:δd=>δd)

show(to)

element = getElements(nodes, entities["Ω"], 10)
prescribe!(element,:EA=>EA,:ρA=>ρA)
prescribe!(element,:u=>(x,y,z)->𝑢(x,y))
prescribe!(element,:δu=>(x,y,z)->δ𝑢(x,y))
prescribe!(element,:∂u∂x=>(x,y,z)->𝑢ₓ(x,y))
prescribe!(element,:∂u∂t=>(x,y,z)->𝑢̇(x,y))
prescribe!(element,:∂δu∂x=>(x,y,z)->δ𝑢ₓ(x,y))
prescribe!(element,:∂δu∂t=>(x,y,z)->δ𝑢̇(x,y))
set∇𝝭!(element)
Hₑu, Hₑδu = Hₑ(element)
println(Hₑu)
println(Hₑδu)
println(log10(Hₑu))
# ──────────────────────────────────────────────────────────
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y/b
    points[3,i] = 0.0
end

# VTK_TRIANGLE = 5,
# VTK_TRIANGLE_STRIP = 6
# VTK_QUAD = 9

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[x.𝐼 for x in elm.𝓒]) for elm in element]
# cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[x.𝐼 for x in elm.𝓒]) for elm in element]
# vtk_grid("./vtk/square_$n.vtu",points,cells) do vtk
# vtk_grid("./vtk/square_$nˣ"*"_$nʸ.vtu",points,cells) do vtk
# vtk_grid("./vtk/square_constraint_all_$nˣ"*"_$nʸ.vtu",points,cells) do vtk
vtk_grid("./vtk/$filename.vtu",points,cells) do vtk
    vtk["d"] = [node.d for node in nodes]
    vtk["δd"] = [node.δd for node in nodes]
    vtk["𝑢"] = [𝑢(node.x,node.y) for node in nodes]
    vtk["δ𝑢"] = [δ𝑢(node.x,node.y) for node in nodes]
end
# ──────────────────────────────────────────────────────────

# end
gmsh.finalize()