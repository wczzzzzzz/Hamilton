
using ApproxOperator, LinearAlgebra
include("import_hmd.jl")
ndiv= 4
elements,nodes = import_hmd_bar("./msh/bar/L=4/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
import ApproxOperator.Hamilton: ∫qmpdΩ, ∫qkpdΩ
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ

m = zeros(nₚ,nₚ)
k = zeros(nₚ,nₚ)
prescribe!(elements["Ω"],:EA=>(x,y,z)->1.0)
prescribe!(elements["Ω"],:ρA=>(x,y,z)->1.0)
ａ = ∫qmpdΩ=>elements["Ω"]
b = ∫qkpdΩ=>elements["Ω"]
ａ(m)
b(k)
println("M has zeros? ", all(m .== 0))
println("K has zeros? ", all(k .== 0))
println("M sum: ", sum(m))
println("K sum: ", sum(k))
println(diag(m))
println(diag(k))

