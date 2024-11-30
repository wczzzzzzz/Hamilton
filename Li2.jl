using ApproxOperator, JuMP, Ipopt
using LinearAlgebra, CairoMakie, GLMakie, XLSX

function V(r, D, α, r₀)
    return D * (exp(-2α * (r - r₀)) - 2 * exp(-α * (r - r₀)))
end

function grad_V(r, D, α, r₀)
    return -D * α * (2 * exp(-2α * (r - r₀)) - 2 * exp(-α * (r - r₀)))
end

model = Model(Ipopt.Optimizer)

@variable(model, r₁)
@variable(model, r₂)

D = 5.0 
α = 1.0  
r₀ = 2.0  
r₁_fixed = 1.0 
r₂_fixed = 3.0  

@objective(model, Min, V(r₁, D, α, r₀) + V(r₂, D, α, r₀))

@constraint(model, r₁ == r₁_fixed)
@constraint(model, r₂ == r₂_fixed)

optimize!(model)

r₁_value = value(r₁)
r₂_value = value(r₂)

println("Optimal r₁: ", r₁_value)
println("Optimal r₂: ", r₂_value)

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1, 1], xlabel = "r₁", ylabel = "r₂", title = "Optimal Positions")
CairoMakie.scatter!(ax, [r₁_value], [r₂_value], color = :red)
fig

# XLSX.writetable("results.xlsx", ["r₁" => [r₁_value], "r₂" => [r₂_value]])

GLMakie.meshscatter([r₁_value], [r₂_value], [0.0], marker = 'o', color = :red)