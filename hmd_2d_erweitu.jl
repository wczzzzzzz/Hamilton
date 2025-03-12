using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 20
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/Non-uniform_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/square/square_"*string(ndiv)*".msh")
elements,nodes = import_hmd_Tri3("./msh/test_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Quad("./msh/test_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ωᵍ"])

ρA = 1e0
EA = 1.0
α = 1e15
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
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->-𝑇(y))
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))

𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₄"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
# 𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]
# 𝑎ᵞ = ∫∫∇v∇udxdy=>elements["Ω"][[146,82,59,175,165,71,134,147].-56]

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

dt =(k+kᵅ)\(f+fᵅ)
# dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
δd = dt[nₚ+1:end]

push!(nodes,:d=>d,:δd=>δd)

fig = Figure()
# ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])
ax1 = Axis(fig[1,1])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)

xs = [i for (i, node) in enumerate(nodes) if node.x == 0]
ys = [nodes[i].y for i in xs]
ds = [nodes[i].d for i in xs]
# xs[i] = node.x
# ys[i] = node.y
# ds[i] = node.d

face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
# meshscatter!(ax1,ys,ds,color=ds,markersize = 0.05)
# meshscatter!(ax1, ys, ds, color=ds, markersize = 0.05)
lines!(ax1, ys[[2,3:end...,1]], ds[[2,3:end...,1]], color = ds, linewidth = 2)
fig

# save("./fig/hmd_2d_二维图_Tri3/t=100.png",fig)
# save("./fig/hmd_2d_二维图_Quad/t=100.png",fig)
# save("./fig/锁三边x=20/二维图/t=16.png",fig)


    