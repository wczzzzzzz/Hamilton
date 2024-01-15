
import Gmsh: gmsh

function import_roof_gauss(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 5
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], Element{:Seg2}, integrationorder)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"], Element{:Poi1}, integrationorder)
    push!(elements["Ω"], :𝝭=>:𝑠,:∂𝝭∂x=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠)

    # gmsh.finalize()
    return elements, nodes
end