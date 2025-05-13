
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

    # gmsh.finalize()
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
    elements["Ωᵗ"], nodes_t, edges = Tri3toTriHermite(elements["Ω"], nodes)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    elements["Γ₁ᵗ"] = Seg2toSegHermite(elements["Γ₁"], nodes, edges)
    elements["Γ₂ᵗ"] = Seg2toSegHermite(elements["Γ₂"], nodes, edges)
    elements["Γ₃ᵗ"] = Seg2toSegHermite(elements["Γ₃"], nodes, edges)
    elements["Γ₄ᵗ"] = Seg2toSegHermite(elements["Γ₄"], nodes, edges)
    push!(elements["Ωᵗ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ₁ᵗ"], :𝝭)
    push!(elements["Γ₂ᵗ"], :𝝭)
    push!(elements["Γ₃ᵗ"], :𝝭)
    push!(elements["Γ₄ᵗ"], :𝝭)

    # gmsh.finalize()
    return elements, nodes, nodes_t
end


