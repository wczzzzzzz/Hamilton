using ApproxOperator, CairoMakie, LinearAlgebra, XLSX
import ApproxOperator.Hamilton: ∫∫Fvdxdt_2dof, ∫∫q̇mṗqkpdxdt_2dof
import ApproxOperator.Heat: H₁, L₂

include("import_hmd.jl")

ndiv = 64
elements, nodes = import_hmd_bar("./msh/bar/L=8/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵍ"])

M = 2.0      
m_b = 1.0        
g = 9.8         
θ = 30.0     
θ_rad = deg2rad(θ)

M_mat = [M+m_b          m_b*cos(θ_rad);
         m_b*cos(θ_rad) m_b           ]

k_mat = zeros(2, 2)

F_vec = [0.0;
         m_b * g * sin(θ_rad)]


ApproxOperator.Hamilton.eval(:(M_mat = $M_mat))
ApproxOperator.Hamilton.eval(:(k_mat = $k_mat))
ApproxOperator.Hamilton.eval(:(F_vec = $F_vec))

T_end = 2.0

a_exact = M_mat \ F_vec
a_X = a_exact[1]
a_s = a_exact[2]

q_exact_X(t) = 0.5 * a_X * t^2
q_exact_s(t) = 0.5 * a_s * t^2

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->q_exact_s(x))

k = zeros(2nₚ, 2nₚ)
f = zeros(2nₚ)

𝑎 = ∫∫q̇mṗqkpdxdt_2dof => elements["Ω"]
𝑎(k)

b = ∫∫Fvdxdt_2dof => elements["Ω"]
b(f)

v₀ = [0.0, 0.0]
q₀ = [0.0, 0.0]
𝑃₀ = M_mat * v₀
f[1:2] .-= 𝑃₀  

α = 1e12
kᵅ = zeros(2nₚ, 2nₚ)
fᵅ = zeros(2nₚ)
kᵅ[1:2, 1:2] .+= α * I(2)
fᵅ[1:2] .+= α .* q₀

kᵝ = zeros(2nₚ, 2nₚ)
fᵝ = zeros(2nₚ)
kᵝ[(2nₚ-1):2nₚ, (2nₚ-1):2nₚ] .+= α * I(2)

sol = [k + kᵅ  -k;
       -k       kᵝ] \ [fᵅ; -f + fᵝ]

d = sol[1:2nₚ]

d_s = d[2:2:end]

push!(nodes, :d => d_s)

𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
# println("log10(L2) for s: ", 𝐿₂)
println(𝐿₂)

# index = [8, 16, 32, 64, 128, 256]
# XLSX.openxlsx("./excel/楔块滑块.xlsx", mode="rw") do xf
#     Sheet = xf[1] 
#     ind = findfirst(n->n==ndiv, index) + 1
#     Sheet["A"*string(ind)] = log10(8/ndiv)
#     Sheet["C"*string(ind)] = 𝐿₂
# end
