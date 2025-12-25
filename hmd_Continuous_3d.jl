using  ApproxOperator, XLSX, LinearAlgebra, LinearSolve
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie

include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 10
elements,nodes = import_hmd_3d("./msh/3d/"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri6("./msh/Non-uniform/618/Tri6_2.0_"*string(ndiv)*".msh")

nₚ = length(nodes)
nₑ = length(elements["Ω"])

# set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set𝝭!(elements["Γ₅"])
set𝝭!(elements["Γ₆"])
set∇𝝭!(elements["Ωᵍ"])

α = 1e7
ρA = 1.0
EA = 1.0
a = 2.0^0.5
l = 1.0
c = (EA/ρA)^0.5
φ(x,y) = 2*a*π*sin.(π.*x/l)*sin.(π.*y/l)
𝑢(x,y,t) = 2*sin.(π.*x/l)*sin.(π.*y/l)*sin.(π.*a.*t)
# lines!(t, 𝑥, color = :black)
# ∂u∂t(x,t) = (-π.*a.*c/l)*sin.(π.*a.*c*t/l).*sin.(π.*x/l)
# ∂u∂x(x,t) = (π./l)*cos.(π.*a.*c*t/l).*cos.(π.*x/l)

prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₅"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₆"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Ω"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₅"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₆"],:α=>(x,y,z)->α)
prescribe!(elements["Γ₁"],:t=>(x,y,z)->φ(x,y))
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:c=>(x,y,z)->c)

# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)

k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵞ = zeros(nₚ,nₚ)

𝑎 = ∫∫∇q∇pdxdt=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γ₁"]

# 𝑎ᵅ = ∫vgdΓ=>elements["Γ₂"]∪elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₁"]
𝑎ᵅ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₅"]∪elements["Γ₆"]∪elements["Γ₁"]
𝑎ᵝ = ∫vgdΓ=>elements["Γ₃"]∪elements["Γ₄"]∪elements["Γ₅"]∪elements["Γ₆"]∪elements["Γ₂"]

𝑎ᵝ(kᵝ,fᵝ)
𝑎ᵅ(kᵅ,fᵅ)
𝑓(f)
𝑎(k)

dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# dt = (k+kᵅ)\(f+fᵅ)
d = dt[1:nₚ]
# δd = dt[nₚ+1:end]
push!(nodes,:d=>d)
# push!(nodes,:δd=>δd)

# ed = test_domain_error(elements["Ω"])
# e3 = test_boundary_error(elements["Γ₃ₜ"])
# println(ed)
# println(e3)

# 𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
# 𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
# println(𝐿₂)

# index = [4,8,16,32]
# # index = [5,10,20,40]
# XLSX.openxlsx("./excel/hmd_Continuous(2).xlsx", mode="rw") do xf
#     Sheet = xf[7]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     # Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_Continuous/uniform_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
# end

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_Continuous_3d/uniform_"*string(ndiv), points, cells) do vtk
#     vtk["位移"] = d
#     vtk["精确解"] = us
# end
