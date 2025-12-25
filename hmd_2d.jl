using  ApproxOperator
using  BiRefine
using WriteVTK
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫pudΩ, ∫uudΩ, ∫ppdΩ, stabilization_bar_LSG, stabilization_bar_LSG_Γ, truncation_error, test_boundary_error, test_domain_error
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie, XLSX, LinearAlgebra, LinearSolve

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")
# include("importmsh.jl")

ndiv= 32
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/Tri6_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/618/Tri6_0.5_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/square/square_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/局部加密/C=0.2/Tri3_"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/Tri6/"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform/拉伸压缩/2.1_"*string(ndiv)*".msh");uniform = "nonuniform"
# elements,nodes = import_hmd_Tri3("./msh/square/Tri3反向"*string(ndiv)*".msh");uniform = "uniform"
# elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")

elements,nodes = import_hmd_Tri3("./msh/BiRefine/2d/impact_4_refined_r13.msh");uniform = "uniform"

nₚ = length(nodes)
nₑ = length(elements["Ω"])

# set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
# set∇𝝭!(elements["Γ₃ₜ"])
# set∇𝝭!(elements["Γ₄ₜ"])
set∇𝝭!(elements["Ωᵍ"])

# ρA = 1.0*25.0/100.0
ρA = 1.0
EA = 1.0
α = 1e6
# β = 1e12
c = (EA/ρA)^0.5
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(c*t - x)))/π
    end
end
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (c*t - x))
    end
end
function ∂u∂x(x, t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return -sin(π*(c*t - x))
    end
end
# function ∂²u∂t²(x, t)
#     if x < t - 1 || x > t
#         return 0.0
#     else
#         return π * cos(π * (t - x))
#     end
# end
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Ω"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->-𝑇(y))
# prescribe!(elements["Γ₃ₜ"],:EA=>(x,y,z)->EA)
# prescribe!(elements["Γ₄ₜ"],:EA=>(x,y,z)->EA)
# prescribe!(elements["Γ₃ₜ"],:ρA=>(x,y,z)->ρA)
# prescribe!(elements["Γ₄ₜ"],:ρA=>(x,y,z)->ρA)
# prescribe!(elements["Γ₃ₜ"],:α=>(x,y,z)->α)
# prescribe!(elements["Γ₄ₜ"],:α=>(x,y,z)->α)
prescribe!(elements["Ω"],:c=>(x,y,z)->c)

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)


𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₄"]
# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₁"]∪elements["Γ₂"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₂"]
# 𝑎ᵞ = stabilization_bar_LSG_Γ=>elements["Γ₄ₜ"]∪elements["Γ₃ₜ"]
# 𝑎ᵞ = [
#     stabilization_bar_LSG=>elements["Ω"],
#     # stabilization_bar_LSG_Γ=>elements["Γ₄ₜ"]∪elements["Γ₃ₜ"],
# ]

k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵞ = zeros(nₚ,nₚ)
kᵗ = zeros(nₚ,nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(kᵅ,fᵅ)
𝑎ᵝ(kᵝ,fᵝ)
# 𝑎ᵞ(kᵞ)

# kᵗ = inv(k + kᵅ)
# # kˢ = -k*kᵗ*k' + kᵝ
# kˢ = [k+kᵅ -k;-k kᵝ]
# C = condskeel(kˢ)
# println(C)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt = [k+kᵅ+kᵞ -k-kᵞ;-k-kᵞ kᵝ+kᵞ]\[fᵅ;-f+fᵝ]
# dt =(k+kᵅ)\(f+fᵅ)
# prob = LinearProblem([k+kᵅ+kᵞ -k-kᵞ;-k-kᵞ kᵝ+kᵞ], [fᵅ;-f+fᵝ])
# sol = solve(prob)
# dt = sol.u

d = dt[1:nₚ]
δd = dt[nₚ+1:end]

push!(nodes,:d=>d,:δd=>δd)

# ed = test_domain_error(elements["Ω"])
# e3 = test_boundary_error(elements["Γ₃ₜ"])
# e4 = test_boundary_error(elements["Γ₄ₜ"])
# println(ed)
# println(e3)
# println(e4)

𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
println(𝐻₁,𝐿₂)

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     # Δ = d[i] - 𝑢(x,y)
#         index = [10,20,40,80]
#         XLSX.openxlsx("./excel/square.xlsx", mode="rw") do xf
#         Sheet = xf[3]
#         ind = findfirst(n->n==ndiv,index)+1
#         # Sheet["A"*string(ind)] = x
#         # Sheet["B"*string(ind)] = y
#         # Sheet["C"*string(ind)] = d₁
#         # Sheet["D"*string(ind)] = Δ
#         Sheet["E"*string(ind)] = 𝐿₂
#         Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end

index = [4,8,16,32]
# index = [0.4,0.3,0.2,0.1]
# index = [0,1,2,3]
XLSX.openxlsx("./excel/hmd_BiRefine.xlsx", mode="rw") do xf
    Sheet = xf[3]
    ind = findfirst(n->n==ndiv,index)+1
    Sheet["A"*string(ind)] = log10(4/ndiv)
    # Sheet["A"*string(ind)] = log10(nₚ)
    Sheet["B"*string(ind)] = 𝐻₁
    Sheet["C"*string(ind)] = 𝐿₂
end

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
us = zeros(nₚ)
# qs = zeros(nₚ)
# as = zeros(nₚ)
es = zeros(nₚ)

for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    us[i] = 𝑢(x,y)
    # qs[i] = ∂u∂t(x,y)
    # as[i] = ∂²u∂t²(x,y)
end
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
    δds[i] = node.δd
    es[i] = ds[i] - us[i]
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# # mesh!(ax,xs,ys,zs,face,color=ds)
# # meshscatter!(ax1,xs,ys,us,color=us,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# # meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.06)
fig

# save("./fig/hmd_2d/test_x=20/t=98.png",fig)
# save("./fig/72测试/Tri6_非均布_LSG_32.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri3/三维图/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/均布/t=25.png",fig)
# save("./fig/hmd_2d/局部加密C=0.2/T6_c=0.05.png",fig)
# save("./fig/hmd_2d/Tri3/非均布/n=80.png",fig)


# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d*4
#     # points[3,i] = us[i]*4
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./vtk/hmd_2d/Tri3_d_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
#     # vtk["精确解"] = us
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_2d/error/uniform_Tri3_"*string(ndiv), points, cells) do vtk
#     # vtk["fₓ"] = fₓ
#     # vtk["fₜ"] = fₜ
#     # vtk["fₓₓ"] = fₓₓ
#     # vtk["fₜₜ"] = fₜₜ
#     # vtk["fₓₓ/fₜₜ"] = fₓₓ./fₜₜ
#     vtk["误差"] = es
# end

