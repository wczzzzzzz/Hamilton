using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX, LinearAlgebra
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie


include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 111
# elements,nodes = import_hmd_Tri3("./msh/square_"*string(ndiv)*".msh")
elements,nodes = import_hmd_Tri6("./msh/tri6_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform_Tri6/"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

α = 1e10
ρA = 1e0
EA = 1.0
φ(x) = sin(π*x/l)
𝑢(x,t) = cos.(π.*a.*t/l).*sin.(π.*x/l)

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

k = zeros(nₚ,nₚ)
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

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
zs = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    # zs[i] = 𝑢(x,y)
    ds[i] = node.d
    # δds[i] = node.δd
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.1)
fig


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

# save("./fig/连续解/锁时间末端Tri_6非均布/t=19.png",fig)
# save("./fig/连续解/锁时间末端Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6均布/t=25.png",fig)
# save("./fig/连续解/mix_Tri_6非均布/t=50.png",fig)
