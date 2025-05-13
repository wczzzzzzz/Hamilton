
using ApproxOperator
import ApproxOperator.Test: cc𝝭, cc∇𝝭

include("importmsh.jl")

elements, nodes, nodes_t, = import_hermite("trih.msh")
nₚ = length(nodes)
nₑ = length(elements["Ωᵗ"])
nₗ = length(nodes_t) - nₚ - nₑ
set𝝭!(elements["Ωᵗ"])
set∇𝝭!(elements["Ωᵗ"])
set𝝭!(elements["Γ₁ᵗ"])
set𝝭!(elements["Γ₂ᵗ"])
set𝝭!(elements["Γ₃ᵗ"])
set𝝭!(elements["Γ₄ᵗ"])

# constant
u(x,y,z) = 1.0
∂u∂x(x,y,z) = 0.0
∂u∂y(x,y,z) = 0.0

# linear
# u(x,y,z) = 1+2x+3y
# ∂u∂x(x,y,z) = 2.0
# ∂u∂y(x,y,z) = 3.0

# quadratic
# u(x,y,z) = 1+2x+3y+4x^2+5x*y+6y^2
# ∂u∂x(x,y,z) = 2.0+8x+5y
# ∂u∂y(x,y,z) = 3.0+12y+5x

# cubic
# u(x,y,z) = 1+2x+3y+4x^2+5x*y+6y^2+7x^3+8x^2*y+9x*y^2+10y^3
# ∂u∂x(x,y,z) = 2.0+8x+5y+21x^2+16x*y+9y^2
# ∂u∂y(x,y,z) = 3.0+12y+5x+8x^2+18x*y+30y^2

prescribe!(elements["Ωᵗ"], :u=>u)
prescribe!(elements["Ωᵗ"], :∂u∂x=>∂u∂x)
prescribe!(elements["Ωᵗ"], :∂u∂y=>∂u∂y)
prescribe!(elements["Γ₁ᵗ"],:u=>u)
prescribe!(elements["Γ₂ᵗ"],:u=>u)
prescribe!(elements["Γ₃ᵗ"],:u=>u)
prescribe!(elements["Γ₄ᵗ"],:u=>u)

# Consistency condition
d = zeros(nₚ+nₗ+nₑ)
for i in 1:nₚ
    I = nodes_t[i].𝐼
    x = nodes_t[i].x
    y = nodes_t[i].y
    z = nodes_t[i].z
    d[I] = u(x,y,z)
end
for i in 1:nₗ
    I = nodes_t[nₚ+i].𝐼
    x = nodes_t[nₚ+i].x
    y = nodes_t[nₚ+i].y
    z = nodes_t[nₚ+i].z
    s₁ = nodes_t[nₚ+i].s₁
    s₂ = nodes_t[nₚ+i].s₂
    d[nₚ+i] = ∂u∂x(x,y,z)*s₁ + ∂u∂y(x,y,z)*s₂
end
for i in 1:nₑ
    I = nodes_t[nₚ+nₗ+i].𝐼
    x = nodes_t[nₚ+nₗ+i].x
    y = nodes_t[nₚ+nₗ+i].y
    z = nodes_t[nₚ+nₗ+i].z
    d[nₚ+nₗ+i] = u(x,y,z)
end
push!(nodes_t,:d=>d)

# err = cc𝝭(elements["Ωᵗ"])
err = cc∇𝝭(elements["Ωᵗ"])
# err = cc𝝭(elements["Ωᵗ"][1:1])