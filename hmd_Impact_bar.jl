using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie
using WriteVTK

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 20
elements,nodes = import_hmd_Tri6("./msh/square/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/tri6_x=20/"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ωᵍ"])

ρA = 1.0
EA = 1.0
α = 1e15
L = 1.0
v₀ = 1.0

prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:t=>(x,y,z)->1.0)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₁"]

𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₄"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₄"]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
push!(nodes,:d=>d)

xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]

# vtk_grid("./vtk/碰撞_"*string(ndiv)*"_"*string(nₚ), points, cells) do vtk
#     vtk["碰撞"] = d
# end

σ = zeros(nₑ)
for (j,p) in enumerate(elements["Ω"])
    σ_ = 0.0
    𝑤_ = 0.0
    for ξ in p.𝓖
        B₁ = ξ[:∂𝝭∂x]
        ε = 0.0
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(p.𝓒)
            ε += B₁[i]*xᵢ.d
        end
        σ_ += EA*ε*𝑤
        𝑤_ += 𝑤
    end
    σ[j] = σ_/𝑤_
end
vtk_grid("./vtk/应力_"*string(ndiv)*"_"*string(nₚ), points, cells) do vtk
    vtk["位移"] = d
    vtk["应力"] = σ
end


# xs = [node.x for node in nodes]
# ys = [node.y for node in nodes]
# ds = [node.d for node in nodes]
# fig = Figure()
# ax = Axis(fig[1, 1], xlabel = "x", ylabel = "t")
# scatter!(ax, xs, ys, color = d, markersize = 10)
# fig

# fig = Figure()
# ax1 = Axis3(fig[1,1])

# xs = zeros(nₚ)
# ys = zeros(nₚ)
# zs = zeros(nₚ)
# ds = zeros(nₚ)
# # σs = zeros(nₚ)
# δds = zeros(nₚ)
# for (i,node) in enumerate(nodes)
#     xs[i] = node.x
#     ys[i] = node.y
#     ds[i] = node.d
#     # σs[i] = node.σ
# end
# face = zeros(nₑ,6)
# for (i,elm) in enumerate(elements["Ω"])
#     face[i,:] .= [x.𝐼 for x in elm.𝓒]
# end

# meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# # meshscatter!(ax1,xs,ys,σs,color=ds,markersize = 0.06)
# fig