using  ApproxOperator

using WriteVTK
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫pudΩ, ∫uudΩ, ∫ppdΩ, stabilization_bar_LSG, truncation_error
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie, XLSX

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")
# include("importmsh.jl")

ndiv= 10
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/拉伸压缩/Tri6_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
elements,nodes = import_hmd_Tri3("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/局部加密/C=0.2/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/RefineMesh_0.5/"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/拉伸压缩/2.1_"*string(ndiv)*".msh");uniform = "nonuniform"
# elements,nodes = import_hmd_Tri3("./msh/square/Tri3反向"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Quad("./msh/test_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")

nₚ = length(nodes)
nₑ = length(elements["Ω"])

function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(t - x)))/π
    end
end
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (t - x))
    end
end
fig = Figure()
ax1 = Axis3(fig[1,1])

xs = zeros(nₚ)
ys = zeros(nₚ)
zs = zeros(nₚ)

for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    xs[i] = x
    ys[i] = y
    zs[i] = 𝑢(x,y)
    # zs[i] = ∂u∂t(x,y)
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,zs,face,color=ds)
meshscatter!(ax1,xs,ys,zs,color=zs,markersize = 0.06)
fig

# save("./fig/hmd_2d/精确解/非均布n=40.png",fig)

q = zeros(nₚ)
d = zeros(nₚ)
for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    # d[i] = 𝑢(x,y)
    q[i] = ∂u∂t(x,y)
end

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = d[i]*4
#     # points[3,i] = q[i]*4
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/exact_solution/uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["位移"] = d
# end

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/exact_solution/q_2d_non_uniform_Tri3_"*uniform*"_"*string(ndiv), points, cells) do vtk
#     # vtk["位移"] = d
#     vtk["速度"] = q
# end











# ind = 101
# xs = 0.0:4.0/(ind-1):4.0
# ys = 0.0:4.0/(ind-1):4.0
# zs = zeros(ind,ind)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         zs[i,j] = 𝑢(x,y)
#          XLSX.openxlsx("./excel/exact_solution.xlsx", mode="rw") do xf
#          Sheet = xf[4]
#          ind = findfirst(n->n==ndiv,10)+(i-1)*101+j
#          Sheet["B"*string(ind)] = zs[i,j]
#         end
#     end
# end

# for i in 1:101
# x = xs[i]
# y = ys[i]
#      XLSX.openxlsx("./excel/exact_solution.xlsx", mode="rw") do xf
#     Sheet = xf[4]
#     ind = findfirst(n->n==ndiv,11)+i
#     Sheet["C"*string(ind)] = x
#     Sheet["D"*string(ind)] = y
# end
# end
