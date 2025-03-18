
using ApproxOperator, GLMakie

import Gmsh: gmsh

ndiv = 20
gmsh.initialize()
# gmsh.open("./msh/test_x=20/"*string(ndiv)*".msh")
# gmsh.open("./msh/square/"*string(ndiv)*".msh")
gmsh.open("./msh/b=2/Tri6非均布"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γ¹"] = getElements(nodes,entities["Γ¹"])
elements["Γ²"] = getElements(nodes,entities["Γ²"])
elements["Γ³"] = getElements(nodes,entities["Γ³"])
elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"])
# elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
# elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"])
# elements["∂Ω"] = elements["Γᵍ"]∪elements["Γᵗ"]
# elements["Γᵉ"] = getElements(nodes,entities["Γᵉ"])
elements["∂Ω"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]

# gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.0, aspect = :data, azimuth = -0.5*pi, elevation = 0.5*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
hidedecorations!(ax)

x =  nodes.x
y = nodes.y
z = 0
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 0.5,
    color = :black
)

elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
    # x = [x.x for x in elm.𝓒[[1,2,3,4]]]
    # y = [x.y for x in elm.𝓒[[1,2,3,4]]]

    lines!(x,y,linestyle = :dash, linewidth = 0.7, color = :black)
end

# # boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x =  [x.x for x in elm.𝓒]
    y =  [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 0.2, color = :black)
end

# save("./fig/square_"*string(ndiv)*".png",f)
# save("./fig/三角形节点网格/Tri6非均布_b=2_"*string(ndiv)*".png",f)


f