using  ApproxOperator, JuMP, Ipopt, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie


include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 32
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/RefineMesh_1.0/"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/拉伸压缩/2.5_"*string(ndiv)*".msh");uniform = "nonuniform"
elements,nodes = import_hmd_Tri6("./msh/Non-uniform/RefineMesh_1.0/Tri6_"*string(ndiv)*".msh");uniform = "uniform"
nₚ = length(nodes)
nₑ = length(elements["Ω"])

# set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ωᵍ"])

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
# ∂²u∂t²(x,t) = -(π.*a./l)*(π.*a./l)*cos.(π.*a.*t/l).*sin.(π.*x/l)

prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:t=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->φ(x))

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)

k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₁"]

# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₁"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₂"]∪elements["Γ₄"]∪elements["Γ₁"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₂"]

𝑎ᵝ(kᵝ,fᵝ)
𝑎ᵅ(kᵅ,fᵅ)
𝑓(f)
𝑎(k)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt = (k+kᵅ)\(f+fᵅ)
d = dt[1:nₚ]

# d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
# δd = dt[nₚ+1:end]
push!(nodes,:d=>d)
# push!(nodes,:δd=>δd)

# 𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
# 𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
# println(𝐿₂)

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
zs = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
es = zeros(nₚ)
us = zeros(nₚ)
# qs = zeros(nₚ)
# as = zeros(nₚ)
for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    us[i] = 𝑢(x,y)
    # qs[i] = ∂u∂t(x,y)
    # as[i] = ∂²u∂t²(x,y)
end
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
    # δds[i] = node.δd
    es[i] = ds[i] - us[i]
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax1,xs,ys,us,color=us,markersize = 0.1)
# fig

# save("./fig/617测试/非均布_32.png",fig)

# save("./fig/连续解/锁时间末端Tri_6非均布/t=19.png",fig)
# save("./fig/连续解/锁时间末端Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6非均布/n=41.png",fig)

# index = [8,16,32,64]
# # index = [0.4,0.3,0.2,0.1]
# # index = [0,1,2,3]
# XLSX.openxlsx("./excel/hmd_Continuous.xlsx", mode="rw") do xf
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
# vtk_grid("./vtk/hmd_Continuous/uniform_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
# end

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_Continuous/error_uniform_"*string(ndiv), points, cells) do vtk
#     vtk["误差"] = es
#     # vtk["二阶导"] = as
# end
