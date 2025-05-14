using  ApproxOperator

using WriteVTK
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫pudΩ, ∫uudΩ, ∫ppdΩ, stabilization_bar_LSG, truncation_error
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie, XLSX

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")
# include("importmsh.jl")

ndiv= 20
elements,nodes,nodes_t = import_hermite("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/局部加密/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/RefineMesh_0.5/"*string(ndiv)*".msh");uniform = "uniform"

nₚ = length(nodes)
nₑ = length(elements["Ωᵗ"])
nₗ = length(nodes_t) - nₚ - nₑ

# set∇²𝝭!(elements["Ωᵗ"])
set𝝭!(elements["Ωᵗ"])
set∇𝝭!(elements["Ωᵗ"])
set𝝭!(elements["Γ₁ᵗ"])
set𝝭!(elements["Γ₂ᵗ"])
set𝝭!(elements["Γ₃ᵗ"])
set𝝭!(elements["Γ₄ᵗ"])

# ρA = 1.0*25.0/100.0
ρA = 1.0
EA = 1.0
α = 1e7
c = (EA/ρA)^0.5
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(t - x)))/π
    end
end

prescribe!(elements["Ωᵗ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωᵗ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ᵗ"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₄ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Γ₄ᵗ"],:t=>(x,y,z)->-𝑇(y))
# prescribe!(elements["Ωᵗ"],:c=>(x,y,z)->c)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ωᵗ"]
𝑓 = ∫vtdΓ=>elements["Γ₄ᵗ"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁ᵗ"]∪elements["Γ₂ᵗ"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ᵗ"]

k = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
kˢ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
f = zeros(nₚ+nₗ+nₑ)
kᵅ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
fᵅ = zeros(nₚ+nₗ+nₑ)
kᵝ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
fᵝ = zeros(nₚ+nₗ+nₑ)

𝑎(k)
𝑓(f)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt =(k+kᵅ)\(f+fᵅ)
# dt = [k -k;-k+kᵅ kᵝ]\[zeros(nₚ);-f+fᵝ+fᵅ]
d = dt[1:nₚ]
δd = dt[nₚ+1:end]

push!(nodes,:d=>d,:δd=>δd)

# 𝐿₂ = log10.(L₂(elements["Ωᵍ"]))

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     # Δ = d[i] - 𝑢(x,y)
#         index = [10,20,40,80]
#         XLSX.openxlsx("./excel/square.xlsx", mode="rw") do xf
#         Sheet = xf[3]
#         ind = findfirst(n->n==ndiv,index)+1
#         # Sheet["A"*string(ind)] = x
#         # Sheet["B"*string(ind)] = y
#         # Sheet["C"*string(ind)] = d₁
#         # Sheet["D"*string(ind)] = Δ
#         Sheet["E"*string(ind)] = 𝐿₂
#         Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end

# index = [10,20,40,80]
# index = [0.4,0.3,0.2,0.1]
# # index = [0,1,2,3,4]
# XLSX.openxlsx("./excel/xinsuanzi.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = 𝐿₂
#     # Sheet["B"*string(ind)] = log10(4/ndiv)
#     Sheet["B"*string(ind)] = log10(nₚ)
# end

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    # zs[i] = 𝑢(xs,ys)
    ds[i] = node.d
    # δds[i] = node.δd
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,zs,face,color=ds)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.1)
fig


# save("./fig/hmd_2d/四边形节点/t=100.png",fig)
# save("./fig/hmd_2d/局部加密C=0.2/T6_c=0.05.png",fig)
# save("./fig/hmd_2d/hermite/Tri3/非均布/n=30.png",fig)
# save("./fig/hmd_2d/hermite/Tri3/均布/n=20.png",fig)

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y*4/3
#     points[3,i] = node.d*4
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/nonuniform/非连续解/Tri6_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# xs = [node.x for node in nodes]'
# ys = [node.y*10/5 for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/junbuceshi/"*uniform*"_"*string(ndiv), points, cells) do vtk
#     # vtk["fₓ"] = fₓ
#     # vtk["fₜ"] = fₜ
#     # vtk["fₓₓ"] = fₓₓ
#     # vtk["fₜₜ"] = fₜₜ
#     # vtk["fₓₓ/fₜₜ"] = fₓₓ./fₜₜ
#     vtk["位移"] = d
# end

