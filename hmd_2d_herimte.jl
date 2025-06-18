using  ApproxOperator

using WriteVTK
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫pudΩ, ∫uudΩ, ∫ppdΩ, stabilization_bar_LSG, truncation_error
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie, XLSX
using SparseArrays
using LinearAlgebra

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")
# include("importmsh.jl")

ndiv= 4
# elements,nodes,nodes_t = import_hermite("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/局部加密/C=0.4/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
elements,nodes,nodes_t = import_hermite("./msh/Non-uniform/RefineMesh_1.0/"*string(ndiv)*".msh");uniform = "uniform"

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
# function ∂u∂t(x, t)
#     if x < t - 1 || x > t
#         return 0.0
#     else
#         return sin(π * (t - x))
#     end
# end
# function ∂u∂x(x, t)
#     if x < t - 1
#         return 0.0
#     elseif x > t
#         return 0.0
#     else
#         return -sin(π*(t - x))
#     end
# end
prescribe!(elements["Ωᵗ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωᵗ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Ωᵗ"],:u=>(x,y,z)->𝑢(x,y))
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

# prescribe!(elements["Ωᵗ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵗ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
# prescribe!(elements["Ωᵗ"],:∂u∂z=>(x,y,z)->0.0)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ωᵗ"]
𝑓 = ∫vtdΓ=>elements["Γ₄ᵗ"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁ᵗ"]∪elements["Γ₂ᵗ"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ᵗ"]∪elements["Γ₂ᵗ"]

# k = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# kˢ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# f = zeros(nₚ+nₗ+nₑ)
# kᵅ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# fᵅ = zeros(nₚ+nₗ+nₑ)
# kᵝ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)
# fᵝ = zeros(nₚ+nₗ+nₑ)
# kᵗ = zeros(nₚ+nₗ+nₑ,nₚ+nₗ+nₑ)

k = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
kˢ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
f = zeros(nₚ + nₗ + nₑ)
kᵅ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵅ = zeros(nₚ + nₗ + nₑ)
kᵝ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵝ = zeros(nₚ + nₗ + nₑ)

𝑎(k)
𝑓(f)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

# kᵗ = inv(k + kᵅ)
# # kˢ = -k*kᵗ*k' + kᵝ
# kˢ = [k+kᵅ -k;-k kᵝ]
# C = condskeel(kˢ)
# println(C)

# dt = sparse([k+kᵅ -k;-k kᵝ])\[fᵅ;-f+fᵝ]
dt = ([k+kᵅ -k;-k kᵝ])\[fᵅ;-f+fᵝ]


# d = dt[1:nₚ+nₗ+nₑ]
d = dt[1:nₚ]
# δd = dt[nₚ+nₗ+nₑ+1:end]

push!(nodes,:d=>d,:δd=>δd)

# push!(nodes_t,:d=>d,:δd=>δd)
# 𝐿₂ = log10.(L₂(elements["Ωᵗ"]))
# 𝐻₁, 𝐿₂ = log10.(H₁(elements["Ωᵗ"]))

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
# index = [0.4,0.2,0.1,0.05]
# index = [8,16,32,64]
# XLSX.openxlsx("./excel/hermite.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
# es = zeros(nₚ)
# us = zeros(nₚ)

for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    # us[i] = 𝑢(x,y)
end

for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d

    # δds[i] = node.δd
    # es[i] = ds[i] - us[i]
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# # mesh!(ax,xs,ys,zs,face,color=ds)
# # meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.06)
fig

# save("./fig/hmd_2d/四边形节点/t=100.png",fig)
# save("./fig/hmd_2d/局部加密C=0.2/hermite/c=0.02.png",fig)
# save("./fig/hmd_2d/hermite/Tri3/非均布/c=0.1.png",fig)
# save("./fig/hmd_2d/hermite/Tri3/均布/n=20.png",fig)

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d*4
#     # points[3,i] = es[i]
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_2d_hermite/non-uniform_C=0.2_c=_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
#     # vtk["误差"] = es
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_2d_hermite/error/uniform_"*string(ndiv), points, cells) do vtk
#     # vtk["fₓ"] = fₓ
#     # vtk["fₜ"] = fₜ
#     # vtk["fₓₓ"] = fₓₓ
#     # vtk["fₜₜ"] = fₜₜ
#     # vtk["fₓₓ/fₜₜ"] = fₓₓ./fₜₜ
#     vtk["误差"] = es
# end

