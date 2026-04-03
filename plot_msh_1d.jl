using ApproxOperator, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import Gmsh: gmsh

gmsh.initialize()
# gmsh.open("./msh/bar/L=8/bar_128.msh")
gmsh.open("./msh/bar/L=8/bar_un_256.msh")

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"])
elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"])
elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"])
elements["∂Ω"] = elements["Γᵗ"] ∪ elements["Γᵍ"]

f = Figure()

ax = Axis(f[1, 1], 
    # aspect = DataAspect(), 
    xlabel = "",           
    ylabel = "",           
    xticksvisible = false,  
    yticksvisible = false, 
    xticklabelsvisible = false, 
    yticklabelsvisible = false 
)
hidespines!(ax) 
hidedecorations!(ax)

x_nodes = [node.x for node in nodes]
y_nodes = [node.y for node in nodes]
scatter!(ax, x_nodes, y_nodes,
    marker = :circle,
    markersize = 4.0,
    color = :black
)

for elm in elements["Ω"]
    x = [node.x for node in elm.𝓒]
    y = [node.y for node in elm.𝓒]
    lines!(ax, x, y,
        linewidth = 2.0,
        color = :blue
    )
end

for elm in elements["∂Ω"]
    x = [node.x for node in elm.𝓒]
    y = [node.y for node in elm.𝓒]
    scatter!(ax, x, y,
        marker = :star,
        markersize = 2.0,
        color = :black
    )
end

display(f)
# save("./fig/弹簧小车网格图/un_n=256.png",f)
# save("./fig/弹簧小车网格图/n=128.png",f)
