using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, stabilization_bar_LSG, ∫q∇𝑛pds
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv = 15
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/square/square_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/test_x=20/"*string(ndiv)*".msh")
elements,nodes = import_hmd_Tri6("./msh/tri6_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Quad("./msh/test_x=20/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

# set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set∇𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ω∩Γ₃"])
set∇𝝭!(elements["Ω∩Γ₄"])

ρA = 1e0
EA = 1.0
α = 1e8
αₚ = 1e8
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
prescribe!(elements["Γ₃"],:α=>(x,y,z)->αₚ)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:g=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->-𝑇(y))
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ω∩Γ₃"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω∩Γ₃"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Ω∩Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω∩Γ₄"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω∩Γ₄"],:ρA=>(x,y,z)->ρA)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₄"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]
# 𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₂"]
𝑎ᵖ = ∫q∇𝑛pds=>elements["Ω∩Γ₃"]


β = -1e0
k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
kᵖ = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
fᵖ = zeros(nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)
𝑎ᵖ(kᵖ,fᵖ)

dt = [k+kᵅ -k-kᵖ';-k-kᵖ kᵝ]\[fᵅ;-f+fᵝ]
# dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt =(k+kᵅ)\(f+fᵅ)
d = dt[1:nₚ]
δd = dt[nₚ+1:end]

push!(nodes,:d=>d,:δd=>δd)

𝐿₂ = log10(L₂(elements["Ωᵍ"]))

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     # Δ = d[i] - 𝑢(x,y)
#         index = [20,40,80,100]
#         XLSX.openxlsx("./excel/test_x=20.xlsx", mode="rw") do xf
#         Sheet = xf[1]
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
# XLSX.openxlsx("./excel/square.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = 𝐿₂
#     Sheet["B"*string(ind)] = log10(4/ndiv)
# end

# push!(nodes,:d=>d₁)
# for (i,a) in enumerate(elements["Ω"])
#     node1 = a.𝓒[1]
#     node2 = a.𝓒[2]
#     node3 = a.𝓒[3]
#     x1 = node1.x
#     x2 = node2.x
#     x3 = node3.x
#     y1 = node1.y
#     y2 = node2.y
#     y3 = node3.y
#     d1 = node1.d
#     d2 = node2.d
#     d3 = node3.d
#     Δ1 = d1 - 𝑢(x1,y1)
#     Δ2 = d2 - 𝑢(x2,y2)
#     Δ3 = d3 - 𝑢(x3,y3)
#     XLSX.openxlsx("./excel/hmd_2d_n=10.xlsx", mode="rw") do xf
#         Sheet = xf[1]
#         ind = findfirst(n->n==ndiv,11)+6*(i-1)
#         Sheet["A"*string(ind+1)] = x1
#         Sheet["B"*string(ind+1)] = y1
#         Sheet["C"*string(ind+1)] = d1
#         Sheet["D"*string(ind+1)] = Δ1

#         Sheet["A"*string(ind+2)] = x2
#         Sheet["B"*string(ind+2)] = y2
#         Sheet["C"*string(ind+2)] = d2
#         Sheet["D"*string(ind+2)] = Δ2

#         Sheet["A"*string(ind+3)] = x3
#         Sheet["B"*string(ind+3)] = y3
#         Sheet["C"*string(ind+3)] = d3
#         Sheet["D"*string(ind+3)] = Δ3

#         Sheet["A"*string(ind+4)] = 0.5*(x1+x2)
#         Sheet["B"*string(ind+4)] = 0.5*(y1+y2)
#         Sheet["C"*string(ind+4)] = 0.5*(d1+d2)
#         Sheet["D"*string(ind+4)] = 0.5*(Δ1+Δ2)

#         Sheet["A"*string(ind+5)] = 0.5*(x2+x3)
#         Sheet["B"*string(ind+5)] = 0.5*(y2+y3)
#         Sheet["C"*string(ind+5)] = 0.5*(d2+d3)
#         Sheet["D"*string(ind+5)] = 0.5*(Δ2+Δ3)

#         Sheet["A"*string(ind+6)] = 0.5*(x3+x1)
#         Sheet["B"*string(ind+6)] = 0.5*(y3+y1)
#         Sheet["C"*string(ind+6)] = 0.5*(d3+d1)
#         Sheet["D"*string(ind+6)] = 0.5*(Δ3+Δ1)
#     end
# end

# push!(nodes,:d=>d)
# fig = Figure()
# Axis(fig[1, 1])
# xs = [node.x for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# ys = [node.d for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# lines!(xs,ys, color = :blue)

# fig

# ind = 121
# xs = 0.0:4.0/(ind-1):4.0
# ys = 0.0:4.0/(ind-1):4.0
# zs = zeros(ind,ind)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         zs[i,j] = 𝑢(x,y)
#     end
# end

fig = Figure()
ax1 = Axis3(fig[1,1])
ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    # zs[i] = 𝑢(xs,ys)
    ds[i] = node.d
    δds[i] = node.δd
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.05)
meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.05)
# zlims!(ax1, -1, 1)
zlims!(ax2, -1, 1)
fig

# save("./fig/hmd_2d/test_x=20/t=98.png",fig)
# save("./fig/hmd_2d/四边形节点/t=100.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri3/三维图/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/均布/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/非均布/t=15.png",fig)
# save("./fig/hmd_2d/Tri6/均布/t=25.png",fig)


    