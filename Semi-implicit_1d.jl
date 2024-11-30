using  ApproxOperator, JuMP, Ipopt, CairoMakie

model = Model(Ipopt.Optimizer)

include("import_hmd.jl")

ndiv= 80
elements,nodes = import_hmd_bar("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
Axis(fig[1, 1])
𝑡 = 0.0:0.01:1.0
𝜔 = (k/m)^0.5
𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(𝑡, 𝑥, color = :black)

T = 10
Δt = 0.1
nₜ = Int(T/Δt)

q = zeros(nₜ+1)
q̇ = zeros(nₜ+1)
q̈ = zeros(nₜ+1)

q[1] = q₀
q̇[1] = q̇₀

for n in 1:nₜ
    q̈[n+1] = -k/m * q[n]  
    q̇[n+1] = q̇[n] + Δt * q̈[n+1]  
    q[n+1] = q[n] + Δt * q̇[n+1]  
end

lines!(t, q, color = :blue)


fig

# save("./fig/一维/hmd_1d.png",fig)
