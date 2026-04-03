using  ApproxOperator, CairoMakie, XLSX
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫q̇mṗqkpdxdt
import ApproxOperator.Heat: H₁, L₂

# model = Model(Ipopt.Optimizer)

include("import_hmd.jl")

ndiv= 64
elements,nodes = import_hmd_bar("./msh/bar/L=8/bar_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_bar("./msh/bar/L=8/bar_un_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵍ"])

kᶜ = 50.0
m = 1.0
q̇₀ = 1.0
q₀ = 1.0
prescribe!(elements["Ω"],:m=>(x,y,z)->m)
prescribe!(elements["Ω"],:k=>(x,y,z)->kᶜ)


fig = Figure()
Axis(fig[1, 1])
t = 0.0:0.03125:8.0
𝜔 = (kᶜ/m)^0.5
𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
∂u∂x(t) = 0.0
∂u∂t(t) = q̇₀*cos(𝜔*t) - 𝜔*q₀*sin(𝜔*t)

𝑥 = q₀.*cos.(𝜔.*t) + q̇₀/𝜔.*sin.(𝜔.*t)
lines!(t, 𝑥, color = :black)

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x))
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

𝑎 = ∫∫q̇mṗqkpdxdt=>elements["Ω"]

𝑎(k)

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

α = 1e12
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*q₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
# kᵝ[1,1] += α
kᵝ[2,2] += α

d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = d[1:nₚ]

push!(nodes,:d=>d)
𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
println(𝐿₂)
# 𝐻₁,𝐿₂ = log10.(H₁(elements["Ωᵍ"]))
# println(𝐻₁,𝐿₂)

# lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2]], color = :blue)
# lines!(nodes.x[[1,3,27,18,19,17,26,16,33,15,25,14,32,13,24,12,31,11,23,10,30,9,22,8,29,7,21,6,28,5,20,4,2]], 
# d[[1,3,27,18,19,17,26,16,33,15,25,14,32,13,24,12,31,11,23,10,30,9,22,8,29,7,21,6,28,5,20,4,2]], color = :blue)
# lines!(nodes.x, d, color = :blue)

# e = d[[1,3:end...,2]] - 𝑥
# lines!(𝑡, e[[1,3:end...,2]], color = :red)

# fig

# save("./fig/一维/位移图/64.png",fig)

# index = [32,64,128,256]
# XLSX.openxlsx("./excel/弹簧小车.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(8/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end
