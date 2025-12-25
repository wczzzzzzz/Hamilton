using  ApproxOperator, XLSX, LinearAlgebra, LinearSolve
import ApproxOperator.Hamilton: ∫qmpdΩ, ∫qkpdΩ
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie
import Gmsh: gmsh

include("import_hmd.jl")

ndiv= 20
elements,nodes = import_hmd_bar("./msh/bar/L=4/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

α = 1e7
ρA = 1.0
EA = 1.0

k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
fᵗ = zeros(nₚ)
fᵍ = zeros(nₚ)

prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->α)

𝑎 = ∫qmpdΩ=>elements["Ω"]
b = ∫qkpdΩ=>elements["Ω"]
𝑎ᵅ = ∫vgdΓ=>elements["Γᵍ"]

𝑎(m)
b(k)
𝑎ᵅ(m,fᵍ)

T = 4
Δt = 0.0001
nₜ = Int(T/Δt)
d = zeros(nₚ,nₜ+1)
d̈ₙ₊₁ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

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

for n in 1:nₜ
    fill!(fᵗ,0.0)
    t = (n+1)*Δt
    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->-𝑇(t))
    𝑓 = ∫vtdΓ=>elements["Γᵗ"]
    𝑓(fᵗ)
    d̈ₙ₊₁ .= m\(fᵗ+fᵍ - k*d[:,n])
    ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ₊₁
    d[:,n+1] .= d[:,n] + Δt*ḋₙ₊₁
end

push!(nodes,:d=>d[:,end])
# push!(nodes,:d=>d[:,17])

points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = 2.0
    points[3,i] = node.d
end

cells = [MeshCell(VTKCellTypes.VTK_POLY_LINE,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
vtk_grid("./vtk/semi_2d/"*string(ndiv)*".vtu",points,cells) do vtk
    # vtk["d"] = [node.d for node in nodes]
    vtk["位移"] = d[:,end]
end




# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x,y))
# prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)
# 𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
# println(𝐻₁,𝐿₂)

# ys = 0.0:4.0/(17-1):4.0
# for (i, node) in enumerate(nodes)
#     for (j, t) in enumerate(ys)
#         x = node.x
#         z = d[i,j]
#         Δ = d[i,j] - 𝑢(x,t)
#             # index = [8,16,32,64]
#             XLSX.openxlsx("./excel/Semi-implicit_Euler.xlsx", mode="rw") do xf
#             Sheet = xf[3]
#             ind = findfirst(n->n==ndiv,16)+(i-1)*33+j
#             # ind = findfirst(n->n==ndiv,index)+1
#             Sheet["A"*string(ind)] = x
#             Sheet["B"*string(ind)] = t
#             Sheet["C"*string(ind)] = z
#             Sheet["D"*string(ind)] = Δ
#             # Sheet["E"*string(ind)] = log10(4/ndiv)
#             # Sheet["F"*string(ind)] = 𝐻₁
#             # Sheet["G"*string(ind)] = 𝐿₂
#         end
#     end
# end
