
using  ApproxOperator

using WriteVTK
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫pudΩ, ∫uudΩ, ∫ppdΩ, stabilization_bar_LSG, truncation_error
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie, XLSX

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 20
ndivs= 16
elements,nodes,nodes_s = import_hmd_mix("./msh/square/square_"*string(ndiv)*".msh","./msh/square/square_"*string(ndivs)*".msh",ndivs)
nₚ = length(nodes)
nₜ = length(nodes_s)
nₑ = length(elements["Ω"])

# set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
# set∇𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωₚ"])
set𝝭!(elements["Γ₁ₚ"])
set𝝭!(elements["Γ₂ₚ"])
set𝝭!(elements["Γ₃ₚ"])
set𝝭!(elements["Γ₄ₚ"])

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
function P(x,t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return ρA*sin(π*(t - x))
    end
end
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ωₚ"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Γ₃"],:g=>(x,y,z)->𝑢(x,y))
# prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->-𝑇(y))
prescribe!(elements["Γ₃"],:t=>(x,y,z)->P(x,y))
# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ω"],:c=>(x,y,z)->c)

𝑎ₚᵤ = ∫pudΩ=>(elements["Ωₚ"],elements["Ω"])
𝑎ₚₚ = ∫ppdΩ=>elements["Ωₚ"]
𝑎ᵤᵤ = ∫uudΩ=>elements["Ω"]
𝑓₁ = ∫vtdΓ=>elements["Γ₃"]
# 𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₄"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
# 𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]

kₚᵤ = zeros(nₜ,nₚ)
kₚₚ = zeros(nₜ,nₜ)
kᵤᵤ = zeros(nₚ,nₚ)
k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
f = zeros(nₚ)
f₁ = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₜ,nₜ)
fᵝ = zeros(nₜ)

𝑎ₚᵤ(kₚᵤ)
𝑎ₚₚ(kₚₚ)
𝑎ᵤᵤ(kᵤᵤ)

# 𝑎(k)
𝑓(f)
𝑓₁(f₁)
𝑎ᵅ(kᵅ,fᵅ)
prescribe!(elements["Γ₁ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄ₚ"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁ₚ"],:g=>(x,y,z)->P(x,y))
prescribe!(elements["Γ₂ₚ"],:g=>(x,y,z)->P(x,y))
prescribe!(elements["Γ₃ₚ"],:g=>(x,y,z)->P(x,y))
prescribe!(elements["Γ₄ₚ"],:g=>(x,y,z)->P(x,y))
𝑎ᵝ = ∫vgdΓ=>elements["Γ₁ₚ"]∪elements["Γ₂ₚ"]∪elements["Γ₃ₚ"]∪elements["Γ₄ₚ"]

𝑎ᵝ(kᵝ,fᵝ)

dt = [kᵤᵤ+kᵅ kₚᵤ';kₚᵤ kₚₚ+kᵝ]\[fᵅ;fᵝ]
# dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt =(k+kᵅ)\(f+fᵅ)
# dt = [k -k;-k+kᵅ kᵝ]\[zeros(nₚ);-f+fᵝ+fᵅ]
d = dt[1:nₚ]
δd = dt[nₚ+1:end]

push!(nodes,:d=>d,:δd=>δd)


fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    # zs[i] = 𝑢(xs,ys)
    ds[i] = node.d
    # δds[i] = node.δd
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,zs,face,color=ds)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.1)
fig

# save("./fig/hmd_2d/test_x=20/t=98.png",fig)
# save("./fig/hmd_2d/四边形节点/t=100.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri3/三维图/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/均布/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/非均布/t=15.png",fig)
# save("./fig/hmd_2d/Tri6/均布/t=25.png",fig)
