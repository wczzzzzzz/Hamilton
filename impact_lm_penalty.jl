using ApproxOperator, LinearAlgebra 
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫αqṗdxdt, stabilization_bar_LSG, stabilization_bar_LSG_Γ

using TimerOutputs, WriteVTK
import Gmsh: gmsh

ρA = 1.0
EA = 1.0
α = 1e6
c = (EA/ρA)^0.5
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(c*t - x)))/π
    end
end
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (c*t - x))
    end
end
function ∂u∂x(x, t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return -sin(π*(c*t - x))
    end
end

integration_order = 2
const to = TimerOutput()
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/square_tri3_8.msh")
@timeit to "open msh file" gmsh.open("./msh/impact_4_refined_r13.msh")
# @timeit to "open msh file" gmsh.open("./msh/impact_4.msh")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements_Ω = getElements(nodes, entities["Ω"],integration_order)
    prescribe!(elements_Ω, :EA=>EA, :ρA=>ρA)
    @timeit to "calculate shape functions" set∇𝝭!(elements_Ω)
    𝑎 = ∫∫∇q∇pdxdt=>elements_Ω
    @timeit to "assemble" 𝑎(k)
end

@timeit to "calculate ∫vtdΓ on Γ¹" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ¹"],integration_order)
    prescribe!(elements, :g=>0.0, :α=>α)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑎ᵅ = ∫vgdΓ=>elements
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
end

@timeit to "calculate ∫vtdΓ on Γ²" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ²"],integration_order)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑎ᵅ = ∫vgdΓ=>elements
    prescribe!(elements, :g=>0.0, :α=>α)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    prescribe!(elements, :g=>0.0, :α=>α)
    @timeit to "assemble" 𝑎ᵅ(kᵝ,fᵝ)
end

@timeit to "calculate ∫vtdΓ on Γ³" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ³"],integration_order)
    @timeit to "calculate shape functions" set𝝭!(elements)
    prescribe!(elements, :g=>0.0, :α=>α)
    𝑎ᵅ = ∫vgdΓ=>elements
    @timeit to "assemble" 𝑎ᵅ(kᵝ,fᵝ)
end

@timeit to "calculate ∫vtdΓ on Γ⁴" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ⁴"],integration_order)
    prescribe!(elements, :t=>(x,y,z)->𝑇(y))
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑓 = ∫vtdΓ=>elements
    @timeit to "assemble" 𝑓(f)
end

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
δd = dt[nₚ+1:end]
push!(nodes, :d=>d)
push!(nodes,:δd=>δd)


points = zeros(3, nₚ)
for node in nodes
    I = node.𝐼
    points[1,I] = node.x
    points[2,I] = node.y
    points[3,I] = 0.0
    # points[3,I] = node.d
    # points[3,I] = node.δd
end
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements_Ω]
# cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements_Ω]
vtk_grid("vtk/impact.vtu", points, cells) do vtk
    vtk["d"] = [node.d for node in nodes]
    vtk["δd"] = [node.δd for node in nodes]
end


