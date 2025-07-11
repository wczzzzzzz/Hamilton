using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 32
ndiv_p = 16
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform_"*string(ndiv)*".msh")
elements,nodes,nodes_p = import_hmd_mix("./msh/square/square_"*string(ndiv)*".msh","./msh/square/square_"*string(ndiv_p)*".msh",ndiv_p)
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
# set∇𝝭!(elements["Ωᵇ"])
# set∇𝝭!(elements["Ωᵍ"])

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
prescribe!(elements["Ωₚ"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωₚ"],:ρA=>(x,y,z)->ρA)
# prescribe!(elements["Ωᵇ"],:EA=>(x,y,z)->EA)
# prescribe!(elements["Ωᵇ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄ₚ"],:t=>(x,y,z)->-𝑇(y))
# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))

𝑎ᵘ = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑎ᵘᵖ = ∫∫∇q∇pdxdt=>(elements["Ω"],elements["Ωₚ"])
𝑎ᵖ = ∫∫∇q∇pdxdt=>elements["Ωₚ"]
𝑓ᵖ = ∫vtdΓ=>elements["Γ₄ₚ"]
# 𝑎ᵘᵇ = ∫∫∇q∇pdxdt=>(elements["Ω"],elements["Ωᵇ"])
# 𝑎ᵖᵇ = ∫∫∇q∇pdxdt=>(elements["Ωₚ"],elements["Ωᵇ"])
# 𝑎ᵇᵇ = ∫∫∇q∇pdxdt=>elements["Ωᵇ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ₚ"]

kᵤᵤ = zeros(nᵤ,nᵤ)
kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kᵤₚ = zeros(nᵤ,nₚ)
kₚₚ = zeros(nᵤ,nₚ)
fₚ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵘᵇ = zeros(nᵤ,nₑ)
kᵖᵇ = zeros(nₚ,nₑ)
kᵇᵇ = zeros(nₑ,nₑ)

𝑎ᵘ(kᵤᵤ)
𝑎ᵘᵖ(kᵤₚ)
# 𝑎ₚ(kₚₚ)
# 𝑎ᵘᵇ(kᵘᵇ)
# 𝑎ᵖᵇ(kᵖᵇ)
# 𝑎ᵇᵇ(kᵇᵇ)

𝑓ᵖ(fₚ)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

β = 0e0
k̄ᵤᵤ = β*kᵘᵇ*inv(kᵇᵇ)*kᵘᵇ'
k̄ᵤₚ = β*kᵘᵇ*inv(kᵇᵇ)*kᵖᵇ'
k̄ₚₚ = β*kᵖᵇ*inv(kᵇᵇ)*kᵖᵇ'

# dt = [kᵤᵤ+kᵅ -kᵤₚ;-kᵤₚ' kᵝ]\[fᵅ;-fₚ+fᵝ]
dt = [kᵤᵤ-k̄ᵤᵤ+kᵅ -kᵤₚ+k̄ᵤₚ;-kᵤₚ'+k̄ᵤₚ' kᵝ-k̄ₚₚ]\[fᵅ;-fₚ+fᵝ]
d = dt[1:nᵤ]
δd = dt[nᵤ+1:end]

push!(nodes,:d=>d)
push!(nodes_p,:δd=>δd)

# 𝐿₂ = log10(L₂(elements["Ωᵍ"]))

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     Δ = d[i] - 𝑢(x,y)
#         index = [10,20,40,80]
#         XLSX.openxlsx("./excel/Non-uniform.xlsx", mode="rw") do xf
#         Sheet = xf[4]
#         ind = findfirst(n->n==ndiv,index)+i
#         Sheet["A"*string(ind)] = x
#         Sheet["B"*string(ind)] = y
#         Sheet["C"*string(ind)] = d₁
#         Sheet["D"*string(ind)] = Δ
#         # Sheet["E"*string(ind)] = log10(L₂)
#         # Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end


fig = Figure()
ax1 = Axis3(fig[1,1])
ax2 = Axis3(fig[1,2])
# fig

xs = zeros(nᵤ)
ys = zeros(nᵤ)
ds = zeros(nᵤ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    # zs[i] = 𝑢(xs,ys)
    ds[i] = node.d
end
xp = zeros(nₚ)
yp = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes_p)
    xp[i] = node.x
    yp[i] = node.y
    δds[i] = node.δd
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.1)
meshscatter!(ax2,xp,yp,δds,color=δds,markersize = 0.1)
fig

# save("./fig/均布 Γ₁_g_80.png",fig)

    