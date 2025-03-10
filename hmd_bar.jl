using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie
using SparseArrays

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

u = d
v = dt[nₚ+1:end]

function compute_stress(u, EA, dx)
    n = length(u)
    stress = zeros(nₚ)
    for i in 2:n
        ϵ = (u[i] - u[i - 1]) / dx
        stress[i] = EA * ϵ
    end
    stress[1] = stress[2]
    return stress
end

xs = [node.x for node in nodes]

t = 0.5
stress = compute_stress(u, EA, dx)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "Stress")
lines!(ax, xs, stress, color = :blue, linewidth = 2)
fig

# xs = [node.x for node in nodes]
# ys = [node.y for node in nodes]
# ds = [node.d for node in nodes]
# fig = Figure()
# ax = Axis(fig[1, 1], xlabel = "x", ylabel = "t")
# scatter!(ax, ys, ds, color = d, markersize = 10)
# fig

# save("./fig/碰撞/位移图/Tri6均布/t=20.png",fig)


    




# fig = Figure()
# ax1 = Axis3(fig[1,1])

# xs = zeros(nₚ)
# ys = zeros(nₚ)
# zs = zeros(nₚ)
# ds = zeros(nₚ)
# δds = zeros(nₚ)
# for (i,node) in enumerate(nodes)
#     xs[i] = node.x
#     ys[i] = node.y
#     ds[i] = node.d
# end
# face = zeros(nₑ,6)
# for (i,elm) in enumerate(elements["Ω"])
#     face[i,:] .= [x.𝐼 for x in elm.𝓒]
# end

# meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# fig