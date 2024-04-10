
import Gmsh: gmsh

function import_hmd_bar(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 2
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"], integrationorder)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationorder)
    push!(elements["Ω"], :𝝭=>:𝑠,:∂𝝭∂x=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)

    # gmsh.finalize()
    return elements, nodes
end

import Gmsh: gmsh

function import_hmd_Tri3(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 5
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], Element{:Tri3}, integrationorder)
    elements["Γ₁"] = getElements(nodes, entities["Γ₁"], Element{:Seg2}, integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ₂"], Element{:Seg2}, integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ₃"], Element{:Seg2}, integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ₄"], Element{:Seg2}, integrationorder)
    push!(elements["Ω"], :𝝭=>:𝑠,:∂𝝭∂x=>:𝑠,:∂𝝭∂y=>:𝑠)
    push!(elements["Γ₁"], :𝝭=>:𝑠)
    push!(elements["Γ₂"], :𝝭=>:𝑠)
    push!(elements["Γ₃"], :𝝭=>:𝑠)
    push!(elements["Γ₄"], :𝝭=>:𝑠)

    # gmsh.finalize()
    return elements, nodes
end