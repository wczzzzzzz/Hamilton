using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, H₁

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 32
ndiv_p = 32

elements,nodes,nodes_p = import_hmd_mix_uv("./msh/square/Tri6_"*string(ndiv)*".msh",
"./msh/square/square_"*string(ndiv_p)*".msh",ndiv_p)
# elements,nodes,nodes_p = import_hmd_mix_uv("./msh/Non-uniform/Tri6/"*string(ndiv)*".msh",
# "./msh/Non-uniform/Tri3/"*string(ndiv_p)*".msh",ndiv_p)

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

ρA = 1.0
EA = 1.0
α = 1e6
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
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (t - x))
    end
end
function ∂u∂x(x, t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return -sin(π*(t - x))
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
prescribe!(elements["Γ₂ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂ₚ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄ₚ"],:t=>(x,y,z)->-𝑇(y))

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)

𝑎ᵘ = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑎ᵘᵖ = ∫∫∇q∇pdxdt=>(elements["Ω"],elements["Ωₚ"])
𝑓ᵖ = ∫vtdΓ=>elements["Γ₄ₚ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃ₚ"]∪elements["Γ₂ₚ"]

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
# push!(nodes_p,:δd=>δd)

𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
println(𝐿₂)
println(𝐻₁)

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
    # δds[i] = node.δd
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xp,yp,δds,color=δds,markersize = 0.1)
fig

# save("./fig/616测试/非均布_32.png",fig)

# index = [4,8,16,32]
# XLSX.openxlsx("./excel/hmd_2d_mix_uv.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end
 