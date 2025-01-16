using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

model = Model(Ipopt.Optimizer)

include("import_hmd_test.jl")

ndiv= 10
ndivs= 8
elements,nodes,nodes_s = import_hmd_mix("./msh/Non-uniform_"*string(ndiv)*".msh","./msh/Non-uniform_"*string(ndivs)*".msh")
nₚ = length(nodes)
nₜ = length(nodes_s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
set𝝭!(elements["Γ₅"])
set𝝭!(elements["Γ₇"])
set𝝭!(elements["Γ₈"])

α = 1e13
ρA = 1
EA = 1
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₅"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₇"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₈"],:t=>(x,y,z)->𝑇(y))
prescribe!(elements["Γ₄"],:t=>(x,y,z)->𝑇(y))

kₛ = zeros(nₚ,nₜ)
k = zeros(nₚ,nₚ)
f₁ = zeros(nₚ)
f₂ = zeros(nₜ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₜ,nₜ)
fᵝ = zeros(nₜ)

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# kᵅ = zeros(nₚ,nₚ)
# fᵅ = zeros(nₚ)
# kᵝ = zeros(nₚ,nₚ)
# fᵝ = zeros(nₚ)


ops = [
       Operator{:∫∫q̇mpqkpdx}(:ρA=>ρA,:EA=>EA),
       Operator{:∫∫q̇mΨqkΨdx}(:ρA=>ρA,:EA=>EA),
       Operator{:∫𝑃δudx}(),
       Operator{:∫vtdΓ}(),
       Operator{:∫vgdΓ}(:α=>α),
       Operator{:L₂}(),
]



ops[1](elements["Ω"],k)
ops[2](elements["Ω"],elements["Ωˢ"],kₛ)
ops[3](elements["Γ₁"],f₁)
ops[4](elements["Γ₄"],f₁)
ops[3](elements["Γ₅"],f₂)
ops[4](elements["Γ₈"],f₂)
ops[5](elements["Γ₁"],kᵅ,fᵅ)
ops[5](elements["Γ₂"],kᵅ,fᵅ)
ops[5](elements["Γ₇"],kᵝ,fᵝ)


# ops[2](elements["Γ₁"],f)
# ops[3](elements["Γ₄"],f)
# ops[4](elements["Γ₁"],kᵅ,fᵅ)
# ops[4](elements["Γ₂"],kᵅ,fᵅ)
# ops[4](elements["Γ₃"],kᵝ,fᵝ)

# d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
d = [k+kᵅ kₛ;kₛ' kᵝ]\[f₁+fᵅ;f₂+fᵝ]
d₁ = d[1:nₚ]
# d₁ = d[nₚ+1:2nₚ]
# push!(nodes,:d=>d₁)


α = (EA/ρA)^0.5
function 𝑢(x,t)
    if x < α*(t-1)
        return 2*α/π
    elseif α*t < x
        return 0.0
    else
        α/π*(1-cos(π*(t-x/α)))
    end
end

# set𝝭!(elements["Ωᵍ"])
# set∇𝝭!(elements["Ωᵍ"])
# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->𝑢(x,y))
# L₂ = ops[6](elements["Ωᵍ"])

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     Δ = d[i] - 𝑢(x,y)
#         index = 8
#         XLSX.openxlsx("./excel/mix_formulation.xlsx", mode="rw") do xf
#         Sheet = xf[1]
#         ind = findfirst(n->n==ndivs,index)+i
#         Sheet["A"*string(ind)] = x
#         Sheet["B"*string(ind)] = y
#         Sheet["C"*string(ind)] = d₁
#         Sheet["D"*string(ind)] = Δ
#         # Sheet["E"*string(ind)] = log10(L₂)
#         # Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end



    