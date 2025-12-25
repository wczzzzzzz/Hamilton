using  ApproxOperator, JuMP, Ipopt, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie
using SparseArrays

include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 32
elements,nodes,nodes_t = import_hermite("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/局部加密/C=0.2/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri3/"*string(ndiv)*".msh");uniform = "uniform"

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

α = 1e7
ρA = 1.0
EA = 1.0
a = 1.0
l = 4.0
c = (EA/ρA)^0.5
φ(x) = sin(π*x/l)
𝑢(x,t) = cos.(π.*a.*t/l).*sin.(π.*x/l)
∂u∂t(x,t) = (-π.*a./l)*sin.(π.*a.*t/l).*sin.(π.*x/l)
∂u∂x(x,t) = (π./l)*cos.(π.*a.*t/l).*cos.(π.*x/l)

prescribe!(elements["Γ₄ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ᵗ"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Ωᵗ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωᵗ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁ᵗ"],:t=>(x,y,z)->0.0)
prescribe!(elements["Γ₁ᵗ"],:g=>(x,y,z)->φ(x))
prescribe!(elements["Ωᵗ"],:u=>(x,y,z)->𝑢(x,y))

prescribe!(elements["Ωᵗ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵗ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
prescribe!(elements["Ωᵗ"],:∂u∂z=>(x,y,z)->0.0)

k = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
kˢ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
f = zeros(nₚ+nₗ+nₑ)
kᵅ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
fᵅ = zeros(nₚ+nₗ+nₑ)
kᵝ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
fᵝ = zeros(nₚ+nₗ+nₑ)

# k = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
# kˢ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
# f = zeros(nₚ + nₗ + nₑ)
# kᵅ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
# fᵅ = zeros(nₚ + nₗ + nₑ)
# kᵝ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
# fᵝ = zeros(nₚ + nₗ + nₑ)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ωᵗ"]
𝑓 = ∫vtdΓ=>elements["Γ₁ᵗ"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₂ᵗ"]∪elements["Γ₃ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₁ᵗ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₂ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₁ᵗ"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₂ᵗ"]

𝑎ᵝ(kᵝ,fᵝ)
𝑎ᵅ(kᵅ,fᵅ)
𝑓(f)
𝑎(k)

# dt = sparse([k+kᵅ -k;-k kᵝ])\[fᵅ;-f+fᵝ]
dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt = (k+kᵅ)\(f+fᵅ)

d = dt[1:nₚ+nₗ+nₑ]
# d = dt[1:nₚ]

push!(nodes,:d=>d)
push!(nodes_t,:d=>d)

# 𝐿₂ = log10.(L₂(elements["Ωᵗ"]))
𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵗ"]))
# println(𝐿₂)

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
# ds = zeros(nₚ + nₗ + nₑ)
ds = zeros(nₚ)
es = zeros(nₚ)
us = zeros(nₚ)
for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    us[i] = 𝑢(x,y)
end
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
    es[i] = ds[i] - us[i]
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.1)
fig

# save("./fig/连续解/锁时间末端Tri_6非均布/t=19.png",fig)
# save("./fig/连续解/锁时间末端Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6非均布/n=41.png",fig)

# index = [4,8,16,32]
# index = [0.4,0.3,0.2,0.1]
# index = [0,1,2,3]
# XLSX.openxlsx("./excel/hmd_Continuous_hermite.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_Continuous_hermite/uniform_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = d[1:nₚ] 
#     # vtk["误差"] = es
# end

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_Continuous_hermite/error_uniform_"*string(ndiv), points, cells) do vtk
#     vtk["误差"] = es
# end
