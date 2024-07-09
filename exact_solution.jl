using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX

using GLMakie

model = Model(Ipopt.Optimizer)

include("import_hmd_test.jl")

ndiv= 10
elements,nodes = import_hmd_Tri3("./msh/square_"*string(ndiv)*".msh")
nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

α = 1e9
ρA = 1
EA = 1
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:t=>(x,y,z)->𝑇(y))

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

ops = [
       Operator{:∫∫q̇mpqkpdx}(:ρA=>ρA,:EA=>EA),
       Operator{:∫𝑃δudx}(),
       Operator{:∫vtdΓ}(),
       Operator{:∫vgdΓ}(:α=>α),
]



ops[1](elements["Ω"],k)
ops[2](elements["Γ₁"],f)
ops[3](elements["Γ₄"],f)
ops[4](elements["Γ₁"],kᵅ,fᵅ)
ops[4](elements["Γ₂"],kᵅ,fᵅ)
ops[4](elements["Γ₃"],kᵝ,fᵝ)

# d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
# d₁ = d[1:nₚ]
# # δd = d[nₚ+1:end]
# d = d[1:nₚ]


α = (EA/ρA)^0.5
function 𝑢(x,t)
    if x < α*(t-1)
        return 2*α/π
    elseif α*t < x
        return 0
    else
        α/π*(1-cos(π*(t-x/α)))
    end
end

ind = 101
xs = 0.0:4.0/(ind-1):4.0
ys = 0.0:4.0/(ind-1):4.0
zs = zeros(ind,ind)
for (i,x) in enumerate(xs)
    for (j,y) in enumerate(ys)
        zs[i,j] = 𝑢(x,y)
         XLSX.openxlsx("./excel/exact_solution.xlsx", mode="rw") do xf
         Sheet = xf[4]
         ind = findfirst(n->n==ndiv,10)+(i-1)*101+j
         Sheet["B"*string(ind)] = zs[i,j]
        end
    end
end

# for i in 1:101
# x = xs[i]
# y = ys[i]
#      XLSX.openxlsx("./excel/exact_solution.xlsx", mode="rw") do xf
#     Sheet = xf[4]
#     ind = findfirst(n->n==ndiv,11)+i
#     Sheet["C"*string(ind)] = x
#     Sheet["D"*string(ind)] = y
# end
# end
