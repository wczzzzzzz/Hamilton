
import Gmsh: gmsh

function import_bell(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 5
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵗ"], nodes_t = Tri3toTri18(elements["Ω"], nodes)
    push!(elements["Ωᵗ"], :𝝭,:∂𝝭∂x)

    gmsh.finalize()
    return elements, nodes, nodes_t
end

function import_hermite(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 5
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵗ"], nodes_t = Tri3toTriHermite(elements["Ω"], nodes)
    push!(elements["Ωᵗ"], :𝝭,:∂𝝭∂x)

    gmsh.finalize()
    return elements, nodes, nodes_t
end


