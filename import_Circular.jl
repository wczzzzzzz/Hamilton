
using Tensors, BenchmarkExample
import Gmsh: gmsh
function import_Circular(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    # integrationOrder = 2     # Tri3
    integrationOrder = 3     # Quad4 
    integrationOrder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵉ"] = getElements(nodes, entities["Γᵉ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end
function import_Circular_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
    integrationOrder = 2      # Tri3
    # integrationOrder = 3      # Quad4
    integrationOrder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵉ"] = getElements(nodes, entities["Γᵉ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)

    gmsh.open(filename2)
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    s = 2.5*5/ndivs*ones(length(nodes_s))
    push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 1,γ = 2)

    gmsh.open(filename1)
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], type, integrationOrder, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωˢ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    # gmsh.finalize()
    return elements, nodes, nodes_s
end

function import_Circular_quad_RI(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    integrationOrder = 3
    integrationOrder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵉ"] = getElements(nodes, entities["Γᵉ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)
    gmsh.open(filename2)
    integrationOrder_Ωˢ = 0
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], integrationOrder_Ωˢ)
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes, nodes_s
end

prescribeForSSUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵉ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵉ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Ω"],:q=>(x,y,z)->F)
end

prescribeForClampedUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵉ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵉ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γᵉ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Γᵉ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Ω"],:q=>(x,y,z)->F)
end