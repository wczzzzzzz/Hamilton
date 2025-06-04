using  ApproxOperator, JuMP, Ipopt, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie
using SparseArrays

include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 32
elements,nodes,nodes_t = import_hermite("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/局部加密/C=0.2/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/RefineMesh_1.0/"*string(ndiv)*".msh");uniform = "uniform"

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

α = 1e14
ρA = 1.0
EA = 1.0
a = 1.0
l = 4.0
c = (EA/ρA)^0.5
φ(x) = sin(π*x/l)
𝑢(x,t) = cos.(π.*a.*t/l).*sin.(π.*x/l)
𝑇(x) = EA*cos(π*a/l)*sin(π*x/l)

prescribe!(elements["Γ₄ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂ᵗ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ᵗ"],:𝑃=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃ᵗ"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ωᵗ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωᵗ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄ᵗ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁ᵗ"],:t=>(x,y,z)->𝑇(x))
prescribe!(elements["Γ₁ᵗ"],:g=>(x,y,z)->φ(x))
# prescribe!(elements["Ωᵗ"],:u=>(x,y,z)->𝑢(x,y))

# k = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# kˢ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# f = zeros(nₚ+nₗ+nₑ)
# kᵅ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# fᵅ = zeros(nₚ+nₗ+nₑ)
# kᵝ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# fᵝ = zeros(nₚ+nₗ+nₑ)

k = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
kˢ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
f = zeros(nₚ + nₗ + nₑ)
kᵅ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵅ = zeros(nₚ + nₗ + nₑ)
kᵝ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵝ = zeros(nₚ + nₗ + nₑ)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ωᵗ"]
𝑓 = ∫vtdΓ=>elements["Γ₁ᵗ"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₂ᵗ"]∪elements["Γ₃ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₁ᵗ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₂ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₁ᵗ"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ᵗ"]∪elements["Γ₄ᵗ"]∪elements["Γ₂ᵗ"]

𝑎ᵝ(kᵝ,fᵝ)
𝑎ᵅ(kᵅ,fᵅ)
𝑓(f)
𝑎(k)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt = (k+kᵅ)\(f+fᵅ)

d = dt[1:nₚ+nₗ+nₑ]
# d = dt[1:nₚ]

push!(nodes,:d=>d)

# push!(nodes_t,:d=>d)
# 𝐿₂ = log10.(L₂(elements["Ωᵗ"]))

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ + nₗ + nₑ)
ys = zeros(nₚ + nₗ + nₑ)
ds = zeros(nₚ + nₗ + nₑ)

for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.1)
fig

# save("./fig/连续解/锁时间末端Tri_6非均布/t=19.png",fig)
# save("./fig/连续解/锁时间末端Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6非均布/n=41.png",fig)

# index = [4,8,16,32]
# # index = [0.4,0.3,0.2,0.1]
# # index = [0,1,2,3]
# XLSX.openxlsx("./excel/hmd_Continuous_hermite.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     Sheet["B"*string(ind)] = 𝐿₂
# end

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y*10/25
#     points[3,i] = node.d
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/nonuniform/连续解/Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
# end