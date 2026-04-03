using ApproxOperator, CairoMakie, XLSX
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫q̇mṗqkpdxdt
import ApproxOperator.Heat: H₁, L₂

include("import_hmd.jl")

ndiv = 256
elements, nodes = import_hmd_bar("./msh/bar/L=8/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵍ"])

kᶜ = 50.0
m = 1.0
q̇₀ = 0.0
q₀ = 0.0
F₀ = 25.0     
Ω_freq = 5.0   

𝜔 = sqrt(kᶜ/m)
U₀ = F₀ / (kᶜ - m * Ω_freq^2)
A = q₀
B = (q̇₀ - U₀ * Ω_freq) / 𝜔


prescribe!(elements["Ω"],:m=>(x,y,z)->m)
prescribe!(elements["Ω"],:k=>(x,y,z)->kᶜ)
prescribe!(elements["Ω"],:F=>(x,y,z)->F₀*sin(Ω_freq*x)) 


𝑢(t) = A*cos(𝜔*t) + B*sin(𝜔*t) + U₀*sin(Ω_freq*t)
# ∂u∂x(t) = 0.0
# ∂u∂t(t) = -A*𝜔*sin(𝜔*t) + B*𝜔*cos(𝜔*t) + U₀*Ω_freq*cos(Ω_freq*t)

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x))
# prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x))
# prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂t(x))
# prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)


k = zeros(nₚ,nₚ)
f = zeros(nₚ)

𝑎 = ∫∫q̇mṗqkpdxdt=>elements["Ω"]
𝑎(k, f)

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

α = 1e12
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*q₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵝ[2,2] += α

d = [k+kᵅ -k; -k kᵝ] \ [fᵅ; -f+fᵝ]
d = d[1:nₚ]
push!(nodes, :d=>d)


𝐿₂ = log10.(L₂(elements["Ωᵍ"]))
println(𝐿₂)

fig = Figure()
ax = Axis(fig[1, 1], title="Forced Harmonic Oscillator (Mesh)", xlabel="Time (s)", ylabel="Displacement")

t = 0.0:0.03125:8.0
𝑥_exact = 𝑢.(t)
lines!(ax, t, 𝑥_exact, color = :black, linewidth=2, label="Exact Solution")

nodes_sorted_idx = [1; 3:lastindex(nodes.x); 2]
lines!(ax, nodes.x[nodes_sorted_idx], d[nodes_sorted_idx], color = :blue, label="Mesh Numerical")

axislegend(ax, position=:rt)
display(fig)

index = [32, 64, 128, 256]
XLSX.openxlsx("./excel/简谐小车.xlsx", mode="rw") do xf
    Sheet = xf[1] 
    ind = findfirst(n->n==ndiv, index) + 1
    Sheet["A"*string(ind)] = log10(8/ndiv)
    # Sheet["B"*string(ind)] = 𝐻₁
    Sheet["C"*string(ind)] = 𝐿₂
end
