
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
    push!(elements["Ω"], :𝝭,:∂𝝭∂x)
    push!(elements["Γᵍ"], :𝝭)
    push!(elements["Γᵗ"], :𝝭)

    gmsh.finalize()
    return elements, nodes
end

import Gmsh: gmsh

function import_hmd_Tri3(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 2
    integrationorder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationorder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    push!(elements["Ω"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ₁"], :𝝭)
    push!(elements["Γ₂"], :𝝭)
    push!(elements["Γ₃"], :𝝭)
    push!(elements["Γ₄"], :𝝭)

    type = PiecewiseParametric{:Bubble,:Tri3}
    # type = PiecewiseParametric{:Bubble,:Quad}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"], type, integrationorder)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)


    # elements["Ω∩Γ₃"] = getBoundaryGradientElement(elements["Γ₃"],elements["Ω"])

    # gmsh.finalize()
    return elements, nodes
end

import Gmsh: gmsh

function import_hmd_Tri6(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 4
    integrationorder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationorder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    push!(elements["Ω"], :𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂x∂y)
    push!(elements["Ωᵍ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂x∂y)
    push!(elements["Γ₁"], :𝝭)
    push!(elements["Γ₂"], :𝝭)
    push!(elements["Γ₃"], :𝝭)
    push!(elements["Γ₄"], :𝝭)

    # gmsh.finalize()
    return elements, nodes
end

import Gmsh: gmsh

function import_hmd_Quad(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 3
    integrationorder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationorder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    push!(elements["Ω"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ₁"], :𝝭)
    push!(elements["Γ₂"], :𝝭)
    push!(elements["Γ₃"], :𝝭)
    push!(elements["Γ₄"], :𝝭)

    type = PiecewiseParametric{:Bubble,:Tri3}
    # type = PiecewiseParametric{:Bubble,:Quad}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"], type, integrationorder)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)


    # elements["Ω∩Γ₃"] = getBoundaryGradientElement(elements["Γ₃"],elements["Ω"])

    # gmsh.finalize()
    return elements, nodes
end

# function getBoundaryGradientElement(as::Vector{T},bs::Vector{S}) where {T,S}
#     elms = S[]
#     ξ = Float64[]
#     η = Float64[]
#     data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
#     data_b = getfield(bs[1].𝓖[1],:data)
#     data[:x] = data_b[:x]
#     data[:y] = data_b[:y]
#     data[:z] = data_b[:z]
#     data[:ξ] = (2,ξ)
#     data[:η] = (2,η)
#     data[:w] = data_b[:w]
#     data[:𝑤] = data_b[:𝑤]
#     data[:n₁] = data_b[:n₁]
#     data[:n₂] = data_b[:n₂]
#     data[:𝐽] = data_b[:𝐽]
#     G = 0
#     C = 0
#     s = 0
#     for a in as
#         indices_a = [xᵢ.𝐼 for xᵢ in a.𝓒]
#         for b in bs
#             indices_b = [xᵢ.𝐼 for xᵢ in b.𝓒]
#             if indices_a ⊂ indices_b
#                 indices_turn = indexin(a,b)
#                 if isa(b,Element{:Tri3})
#                     C += 1
#                     𝓒 = b.𝓒
#                     𝓖 = 𝑿ₛ[]
#                     for (g,xg) in enumerate(a.𝓖)
#                         if indices_turn == [1,2]
#                             push!(ξ,0.5*(1+xg.ξ))
#                             push!(η,0.0)
#                         elseif indices_turn == [2,3]
#                             push!(ξ,0.5*(1-xg.ξ))
#                             push!(η,0.5*(1+xg.ξ))
#                         else
#                             push!(ξ,0.0)
#                             push!(η,0.5*(1-xg.ξ))
#                         end
#                         G += 1
#                         push!(𝓖,typeof(xg)(𝑔=g,𝐺=G,𝐶=C,𝑠=S),data)
#                         s += 3
#                     end
#                     push!(elms,Element{:Tri3}(𝓒,𝓖))
#                 end
#             end
#         end
#     end
#     return elms
# end

function import_hmd_mix(filename1::String,filename2::String,n::Int)
    gmsh.initialize()
    
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 5,γ = 3)
    s = 1.5*4/n*ones(length(nodes_s))
    push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)

    gmsh.open(filename1)
    integrationorder = 8
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    elements["Ωₚ"] = getElements(nodes_s, entities["Ω"], type, integrationorder, sp)
    elements["Γ₁ₚ"] = getElements(nodes_s, entities["Γ¹"], type, integrationorder, sp)
    elements["Γ₂ₚ"] = getElements(nodes_s, entities["Γ²"], type, integrationorder, sp)
    elements["Γ₃ₚ"] = getElements(nodes_s, entities["Γ³"], type, integrationorder, sp)
    elements["Γ₄ₚ"] = getElements(nodes_s, entities["Γ⁴"], type, integrationorder, sp)

    nₘ=21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)

    push!(elements["Ω"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ₁"], :𝝭)
    push!(elements["Γ₂"], :𝝭)
    push!(elements["Γ₃"], :𝝭)
    push!(elements["Γ₄"], :𝝭)
    push!(elements["Ωₚ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ₁ₚ"], :𝝭)
    push!(elements["Γ₂ₚ"], :𝝭)
    push!(elements["Γ₃ₚ"], :𝝭)
    push!(elements["Γ₄ₚ"], :𝝭)
    push!(elements["Ωₚ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ₁ₚ"], :𝗠=>𝗠)
    push!(elements["Γ₂ₚ"], :𝗠=>𝗠)
    push!(elements["Γ₃ₚ"], :𝗠=>𝗠)
    push!(elements["Γ₄ₚ"], :𝗠=>𝗠)

    # type = PiecewiseParametric{:Bubble,:Tri3}
    # type = PiecewiseParametric{:Bubble,:Quad}
    # elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"], type, integrationorder)
    # push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)

    # gmsh.finalize()
    return elements, nodes, nodes_s
end

function import_hermite(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationorder = 5
    # integrationorder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationorder)
    elements["Ωᵗ"], nodes_t = Tri3toTriHermite(elements["Ω"], nodes)
    # elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationorder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], integrationorder)
    # elements["Γ₁ᵗ"], nodes_t = Tri3toTriHermite(elements["Γ₁"], nodes)
    # elements["Γ₂ᵗ"], nodes_t = Tri3toTriHermite(elements["Γ₂"], nodes)
    # elements["Γ₃ᵗ"], nodes_t = Tri3toTriHermite(elements["Γ₃"], nodes)
    # elements["Γ₄ᵗ"], nodes_t = Tri3toTriHermite(elements["Γ₄"], nodes)
    push!(elements["Ωᵗ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    # push!(elements["Ωᵍ"], :𝝭,:∂𝝭∂x,:∂𝝭∂y)
    # push!(elements["Γ₁ᵗ"], :𝝭)
    # push!(elements["Γ₂ᵗ"], :𝝭)
    # push!(elements["Γ₃ᵗ"], :𝝭)
    # push!(elements["Γ₄ᵗ"], :𝝭)
    push!(elements["Γ₁"], :𝝭)
    push!(elements["Γ₂"], :𝝭)
    push!(elements["Γ₃"], :𝝭)
    push!(elements["Γ₄"], :𝝭)
    
    # gmsh.finalize()
    return elements, nodes, nodes_t
end
