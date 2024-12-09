using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX, LinearAlgebra

using GLMakie

model = Model(Ipopt.Optimizer)

include("import_hmd.jl")
# include("import_hmd_test.jl")

ndiv= 40
# elements,nodes = import_hmd_Tri3("./msh/square_"*string(ndiv)*".msh")
elements,nodes = import_hmd_Tri3("./msh/Non-uniform_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

α = 1e9
ρA = 1
EA = 1
l = 4
a = 1
q̇ = 1.0
φ(x) = sin(π*x/l)
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->0.0)
# prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->φ(x))
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:t=>(x,y,z)->q̇)  
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)


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
# ops[3](elements["Γ₃"],fᵅ)
ops[4](elements["Γ₁"],kᵅ,fᵅ)
ops[4](elements["Γ₂"],kᵅ,fᵅ)
# ops[4](elements["Γ₃"],kᵝ,fᵝ)
ops[4](elements["Γ₄"],kᵅ,fᵅ)

# d = k\f
d = (k+kᵅ)\(f+fᵅ)
# d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
# d₁ = d[1:nₚ]
# # δd = d[nₚ+1:end]
d = d[1:nₚ]
push!(nodes,:d=>d)

𝑢(x,t) = cos.(π.*a.*t/l).*sin.(π.*x/l)


fig = Figure()
ax = Axis3(fig[1,1])

xs = zeros(nₚ)
ys = zeros(nₚ)
zs = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    xs[i] = node.x
    ys[i] = node.y
    zs[i] = 𝑢(x,y)
    ds[i] = node.d
    # δds[i] = node.δd
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax,xs,ys,ds,color=ds,markersize = 0.1)
# meshscatter!(ax,xs,ys,δds,color=δds,markersize = 0.1)
fig


# ind = 101
# xs = 0.0:4.0/(ind-1):4.0
# ys = 0.0:4.0/(ind-1):4.0
# zs = zeros(ind,ind)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         zs[i,j] = 𝑢(x,y)
#          XLSX.openxlsx("./excel/exact_solution.xlsx", mode="rw") do xf
#          Sheet = xf[4]
#          ind = findfirst(n->n==ndiv,10)+(i-1)*101+j
#          Sheet["B"*string(ind)] = zs[i,j]
#         end
#     end
# end

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

# save("./fig/连续解/非均布n=80.png",fig)
