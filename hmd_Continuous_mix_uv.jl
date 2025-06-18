using  ApproxOperator, JuMP, Ipopt, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, H₁

using GLMakie


include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 32
ndiv_p = 32

# elements,nodes,nodes_p = import_hmd_mix_uv("./msh/square/Tri6_"*string(ndiv)*".msh",
# "./msh/square/square_"*string(ndiv_p)*".msh",ndiv_p)
elements,nodes,nodes_p = import_hmd_mix_uv("./msh/Non-uniform/RefineMesh_1.0/Tri6_"*string(ndiv)*".msh",
"./msh/Non-uniform/RefineMesh_1.0/"*string(ndiv_p)*".msh",ndiv_p)

nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑ = length(elements["Ω"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ωₚ"])
set𝝭!(elements["Γ₁ₚ"])
set𝝭!(elements["Γ₂ₚ"])
set𝝭!(elements["Γ₃ₚ"])
set𝝭!(elements["Γ₄ₚ"])

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
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₁"],:t=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->φ(x))

prescribe!(elements["Ωₚ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωₚ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₄ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₁ₚ"],:t=>(x,y,z)->0.0)

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
 prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)

kᵤᵤ = zeros(nᵤ,nᵤ)
kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kᵤₚ = zeros(nᵤ,nₚ)
fₚ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

𝑎ᵘ = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑎ᵘᵖ = ∫∫∇q∇pdxdt=>(elements["Ω"],elements["Ωₚ"])
𝑓ᵖ = ∫vtdΓ=>elements["Γ₁ₚ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₂"]∪elements["Γ₄"]∪elements["Γ₁"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ₚ"]∪elements["Γ₄ₚ"]∪elements["Γ₂ₚ"]

𝑎ᵘ(kᵤᵤ)
𝑎ᵘᵖ(kᵤₚ)
𝑓ᵖ(fₚ)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

dt = [kᵤᵤ+kᵅ -kᵤₚ;-kᵤₚ' kᵝ]\[fᵅ;-fₚ+fᵝ]
d = dt[1:nᵤ]
push!(nodes,:d=>d)
δd = dt[nᵤ+1:end]
# push!(nodes_p,:δd=>δd)

# 𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))


fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nᵤ)
ys = zeros(nᵤ)
ds = zeros(nᵤ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
end
# xp = zeros(nₚ)
# yp = zeros(nₚ)
# δds = zeros(nₚ)
# for (i,node) in enumerate(nodes_p)
#     xp[i] = node.x
#     yp[i] = node.y
#     δds[i] = node.δd
# end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xp,yp,δds,color=δds,markersize = 0.1)
# fig

# save("./fig/617测试/32.png",fig)

# index = [4,8,16,32]
# XLSX.openxlsx("./excel/hmd_Continuous_mix_uv.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end