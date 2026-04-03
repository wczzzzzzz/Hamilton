using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Hamilton: ∫∫∇q∇pdΩdt, Hₑ
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ

using TimerOutputs, WriteVTK
# using LinearSolve
using LinearAlgebra
import Gmsh: gmsh

# c = 1e-0
c = 1.0
# 𝑢(x,y,z) = x + 2y + 3z
𝑢(x,y,t) = sin((x+y)/2^0.5-c*t)
𝑢̇(x,y,t) = -c*cos((x+y)/2^0.5-c*t)
𝑢₁(x,y,t) = cos((x+y)/2^0.5-c*t)/2^0.5
𝑢₂(x,y,t) = cos((x+y)/2^0.5-c*t)/2^0.5

δ𝑢(x,y,t) = 0.0
δ𝑢̇(x,y,t) = 0.0
δ𝑢₁(x,y,t) = 0.0
δ𝑢₂(x,y,t) = 0.0

const to = TimerOutput()
gmsh.initialize()

# filename = "cube_hex_20"
filename = "cube_10"
@timeit to "open msh file" gmsh.open("msh/$filename.msh")
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
@timeit to "calculate ∫∫∇q∇pdΩdt" begin
    @timeit to "get elements" element = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(element, :c²=>c^2)
    @timeit to "calculate shape functions" set∇𝝭!(element)
    𝑎 = ∫∫∇q∇pdΩdt=>element
    @timeit to "assemble" 𝑎(k)
end

α = 1e8
@timeit to "calculate boundary for u" begin
    @timeit to "get elements" element_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" element_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" element_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" element_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    @timeit to "get elements" element_5 = getElements(nodes, entities["Γ⁵"],integrationOrder)
    @timeit to "get elements" element_6 = getElements(nodes, entities["Γ⁶"],integrationOrder)
    prescribe!(element_1,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    prescribe!(element_2,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    prescribe!(element_3,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    prescribe!(element_4,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    prescribe!(element_5,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    prescribe!(element_6,:g=>(x,y,z)->𝑢(x,y,z), :α=>α)
    @timeit to "calculate shape functions" set𝝭!(element_1)
    @timeit to "calculate shape functions" set𝝭!(element_2)
    @timeit to "calculate shape functions" set𝝭!(element_3)
    @timeit to "calculate shape functions" set𝝭!(element_4)
    @timeit to "calculate shape functions" set𝝭!(element_5)
    @timeit to "calculate shape functions" set𝝭!(element_6)
    𝑎ᵅ = ∫vgdΓ=>element_1∪element_3∪element_4∪element_5∪element_6
    # 𝑎ᵅ = ∫vgdΓ=>element_1∪element_2∪element_3∪element_4∪element_5∪element_6
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
end
@timeit to "calculate δu" begin
    @timeit to "get elements" element_1 = getElements(nodes, entities["Γ¹"],integrationOrder,normal=true)
    @timeit to "get elements" element_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" element_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" element_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    @timeit to "get elements" element_5 = getElements(nodes, entities["Γ⁵"],integrationOrder)
    @timeit to "get elements" element_6 = getElements(nodes, entities["Γ⁶"],integrationOrder)
    prescribe!(element_1,:t=>(x,y,z,n₁,n₂,n₃)->𝑢̇(x,y,z)*n₃)
    prescribe!(element_1,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    prescribe!(element_2,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    prescribe!(element_3,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    prescribe!(element_4,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    prescribe!(element_5,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    prescribe!(element_6,:g=>(x,y,z)->δ𝑢(x,y,z), :α=>α)
    @timeit to "calculate shape functions" set𝝭!(element_1)
    @timeit to "calculate shape functions" set𝝭!(element_2)
    @timeit to "calculate shape functions" set𝝭!(element_3)
    @timeit to "calculate shape functions" set𝝭!(element_4)
    @timeit to "calculate shape functions" set𝝭!(element_5)
    @timeit to "calculate shape functions" set𝝭!(element_6)
    𝑓 = ∫vtdΓ=>element_1
    𝑎ᵝ = ∫vgdΓ=>element_2∪element_3∪element_4∪element_5∪element_6
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
prescribe!(element,:c²=>c^2)
prescribe!(element,:u=>(x,y,z)->𝑢(x,y,z))
prescribe!(element,:δu=>(x,y,z)->δ𝑢(x,y,z))
prescribe!(element,:∂u∂x=>(x,y,z)->𝑢₁(x,y,z))
prescribe!(element,:∂u∂y=>(x,y,z)->𝑢₂(x,y,z))
prescribe!(element,:∂u∂t=>(x,y,z)->𝑢̇(x,y,z))
prescribe!(element,:∂δu∂x=>(x,y,z)->δ𝑢₁(x,y,z))
prescribe!(element,:∂δu∂y=>(x,y,z)->δ𝑢₂(x,y,z))
prescribe!(element,:∂δu∂t=>(x,y,z)->δ𝑢̇(x,y,z))
set∇𝝭!(element)
Hₑu, Hₑδu = Hₑ(element)
println(Hₑu)
println(Hₑδu)
println(log10(Hₑu))
# ──────────────────────────────────────────────────────────
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = node.z
end

# VTK_TRIANGLE = 5,
# VTK_TRIANGLE_STRIP = 6
# VTK_QUAD = 9
# VTK_TETRA = 10
# VTK_HEXAHEDRON = 12
# VTK_QUADRATIC_TETRA = 24


# cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON,[x.𝐼 for x in elm.𝓒]) for elm in element]
cells = [MeshCell(VTKCellTypes.VTK_TETRA,[x.𝐼 for x in elm.𝓒]) for elm in element]
# cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TETRA,[x.𝐼 for x in elm.𝓒]) for elm in element]
vtk_grid("./vtk/$filename.vtu",points,cells) do vtk
    vtk["d"] = [node.d for node in nodes]
    vtk["δd"] = [node.δd for node in nodes]
    vtk["𝑢"] = [𝑢(node.x,node.y,node.z) for node in nodes]
    vtk["δ𝑢"] = [δ𝑢(node.x,node.y,node.z) for node in nodes]
end
# ──────────────────────────────────────────────────────────

# end
gmsh.finalize()

# check
# d = [𝑢(xᵢ.x,xᵢ.y,xᵢ.z) for xᵢ in nodes]
# err = k*d - f
# println(norm(err))