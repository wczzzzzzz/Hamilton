
using  ApproxOperator, XLSX, LinearAlgebra, LinearSolve
import ApproxOperator.Hamilton: ∫qmpdΩ, ∫qkpdΩ
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁

using GLMakie
import Gmsh: gmsh

include("import_hmd.jl")

ndiv= 20
elements,nodes = import_hmd_bar("./msh/bar/L=4/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

α = 1e9
ρA = 1.0
EA = 1.0
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)

fig = Figure()
Axis(fig[1, 1])

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
Δt = 0.001
nₜ = Int(T/Δt)
d = zeros(nₚ,nₜ+1)
d̈ₙ₊₁ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)

c = (EA/ρA)^0.5
function 𝑢(x,t)
    if x < c*(t-1)
        return 2*c/π
    elseif c*t < x
        return 0
    else
        c/π*(1-cos(π*(t-x/c)))
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