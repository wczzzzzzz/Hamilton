
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

    integrationorder = 2
    integrationorder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], Element{:Tri3}, integrationorder)
    elements["Ω"] = getElements(nodes, entities["Ω"], Element{:Tri3}, integrationorder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], Element{:Tri3}, integrationorder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], Element{:Seg2}, integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], Element{:Seg2}, integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], Element{:Seg2}, integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], Element{:Seg2}, integrationorder)
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


    elements["Ω∩Γ₃"]

    gmsh.finalize()
    return elements, nodes
end

function getBoundaryGradientElement(as::Vector{T},bs::Vector{S}) where {T,S}
    elms = S[]
    x = Float64[]
    y = Float64[]
    z = Float64[]
    ξ = Float64[]
    η = Float64[]
    𝑤 = Float64[]
    n₁ = Float64[]
    n₂ = Float64[]
    𝐽 = Float64[]
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data[:x] = (2,x)
    data[:y] = (2,y)
    data[:z] = (2,z)
    data[:ξ] = (2,ξ)
    data[:η] = (2,η)
    data[:𝑤] = (2,𝑤)
    data[:n₁] = (2,n₁)
    data[:n₂] = (2,n₂)
    for a in as
        indices_a = [xᵢ.𝐼 for xᵢ in a.𝓒]
        for b in bs
            indices_b = [xᵢ.𝐼 for xᵢ in b.𝓒]
            if indices_a ⊂ indices_b
                indices_turn = indexin(a,b)
            end
        end
    end
end

function import_hmd_mix(filename1::String,filename2::String)
    gmsh.initialize()
    
    gmsh.open(filename1)
    integrationorder = 2
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], Element{:Tri3}, integrationorder)
    elements["Γ₁"] = getElements(nodes, entities["Γ¹"], Element{:Seg2}, integrationorder)
    elements["Γ₂"] = getElements(nodes, entities["Γ²"], Element{:Seg2}, integrationorder)
    elements["Γ₃"] = getElements(nodes, entities["Γ³"], Element{:Seg2}, integrationorder)
    elements["Γ₄"] = getElements(nodes, entities["Γ⁴"], Element{:Seg2}, integrationorder)
    
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    s = 2.5*4/ndivs*ones(length(nodes_s))
    push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 1,γ = 2)
    elements["Γ₅"] = getElements(nodes_s, entities["Γ₅"], type, integrationorder, sp)
    elements["Γ₇"] = getElements(nodes_s, entities["Γ₇"], type, integrationorder, sp)
    elements["Γ₈"] = getElements(nodes_s, entities["Γ₈"], type, integrationorder, sp)

    gmsh.open(filename1)
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], type, integrationorder, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))

    push!(elements["Ω"], :𝝭=>:𝑠,:∂𝝭∂x=>:𝑠,:∂𝝭∂y=>:𝑠)
    push!(elements["Ωˢ"], :𝝭=>:𝑠,:∂𝝭∂x=>:𝑠,:∂𝝭∂y=>:𝑠)
    push!(elements["Ωˢ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ₁"], :𝝭=>:𝑠)
    push!(elements["Γ₂"], :𝝭=>:𝑠)
    push!(elements["Γ₃"], :𝝭=>:𝑠)
    push!(elements["Γ₄"], :𝝭=>:𝑠)
    push!(elements["Γ₅"], :𝗠=>𝗠)
    push!(elements["Γ₅"], :𝝭=>:𝑠)
    push!(elements["Γ₇"], :𝝭=>:𝑠)
    push!(elements["Γ₇"], :𝗠=>𝗠)
    push!(elements["Γ₈"], :𝝭=>:𝑠)
    push!(elements["Γ₈"], :𝗠=>𝗠)

    # gmsh.finalize()
    return elements, nodes, nodes_s
end