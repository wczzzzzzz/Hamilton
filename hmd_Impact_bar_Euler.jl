using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using GLMakie
using WriteVTK

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

include("import_hmd.jl")

ndiv= 4
elements,nodes = import_hmd_bar("./msh/bar/L=1/bar_"*string(ndiv)*".msh")

nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

ρA = 1.0
EA = 1.0
α = 1e7
L = 1.0
v₀ = 1.0

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fᵗ = zeros(nₚ)
fᵍ = zeros(nₚ)

prescribe!(elements["Ω"],:EA=>(x,y,z)->EA)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->ρA)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->α)

𝑎 = ∫qmpdΩ=>elements["Ω"]
b = ∫qkpdΩ=>elements["Ω"]
𝑎ᵅ = ∫vgdΓ=>elements["Γᵍ"]

𝑎(m)
b(k)
𝑎ᵅ(m,fᵍ)

T = 2
Δt = 0.25
nₜ = Int(T/Δt)
d = zeros(nₚ,nₜ+1)
d̈ₙ₊₁ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

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

push!(nodes,:d=>d)

xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]

σ = zeros(nₑ)
for (j,p) in enumerate(elements["Ω"])
    σ_ = 0.0
    𝑤_ = 0.0
    for ξ in p.𝓖
        B₁ = ξ[:∂𝝭∂x]
        ε = 0.0
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(p.𝓒)
            ε += B₁[i]*xᵢ.d
        end
        σ_ += EA*ε*𝑤
        𝑤_ += 𝑤
    end
    σ[j] = σ_/𝑤_
end

# vtk_grid("./vtk/impact_bar/Tri6_反向_"*string(ndiv)*"_"*string(nₚ), points, cells) do vtk
#     vtk["位移"] = d
#     vtk["应力"] = σ
# end
