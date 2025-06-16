using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫∇v∇udxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 16
ndiv_p = 16

elements,nodes,nodes_p = import_hmd_mix("./msh/square/Tri6_"*string(ndiv)*".msh","./msh/square/square_"*string(ndiv_p)*".msh",ndiv_p)
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

# set∇𝝭!(elements["Ωᵍ"])

ρA = 1.0
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
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄ₚ"],:t=>(x,y,z)->-𝑇(y))
# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))

𝑎ᵘ = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑎ᵘᵖ = ∫∫∇v∇udxdt=>(elements["Ω"],elements["Ωₚ"])
𝑓ᵖ = ∫vtdΓ=>elements["Γ₄ₚ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ₚ"]

kᵤᵤ = zeros(nᵤ,nᵤ)
kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kᵤₚ = zeros(nᵤ,nₚ)
fₚ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

𝑎ᵘ(kᵤᵤ)
𝑎ᵘᵖ(kᵤₚ)
𝑓ᵖ(fₚ)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)

dt = [kᵤᵤ+kᵅ -kᵤₚ;-kᵤₚ' kᵝ]\[fᵅ;-fₚ+fᵝ]
# dt = [kᵤᵤ-k̄ᵤᵤ+kᵅ -kᵤₚ+k̄ᵤₚ;-kᵤₚ'+k̄ᵤₚ' kᵝ-k̄ₚₚ]\[fᵅ;-fₚ+fᵝ]
d = dt[1:nᵤ]
δd = dt[nᵤ+1:end]

push!(nodes,:d=>d)
push!(nodes_p,:δd=>δd)

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
xp = zeros(nₚ)
yp = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes_p)
    xp[i] = node.x
    yp[i] = node.y
    δds[i] = node.δd
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.1)
# meshscatter!(ax2,xp,yp,δds,color=δds,markersize = 0.1)
fig

# save("./fig/616测试/64.png",fig)

    