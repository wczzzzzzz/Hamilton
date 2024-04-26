using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX

using GLMakie

model = Model(Ipopt.Optimizer)

include("import_hmd_test.jl")

ndiv= 11
elements,nodes = import_hmd_Tri3("./msh/bar_"*string(ndiv)*".msh")
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

d = [k+kᵅ k;k kᵝ]\[f+fᵅ;f+fᵝ]
d₁ = d[1:nₚ]


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


for i in 1:nₚ
    x = nodes.x[i]
    y = nodes.y[i]
    d₁ = d[i]
    Δ = d[i] - 𝑢(x,y)
        XLSX.openxlsx("./excel/hmd_2d_n=10.xlsx", mode="rw") do xf
        Sheet = xf[1]
        ind = findfirst(n->n==ndiv,11)+i
        Sheet["A"*string(ind)] = x
        Sheet["B"*string(ind)] = y
        Sheet["C"*string(ind)] = d₁
        Sheet["D"*string(ind)] = Δ
    end
end


# push!(nodes,:d=>d)
# fig = Figure()
# Axis(fig[1, 1])
# xs = [node.x for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# ys = [node.d for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# lines!(xs,ys, color = :blue)

# fig

# for i = 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     Δ = d[i] - 𝑢(x,y)
#          XLSX.openxlsx("./excel/hmd_2d_error.xlsx", mode="rw") do xf
#             Sheet = xf[1]
#          ind = findfirst(n->n==ndiv,11)+1
#          Sheet["B"*string(ind)] = Δ
#             Sheet = xf[2]
#          ind = findfirst(n->n==ndiv,11)+i
#          Sheet["C"*string(ind)] = x
#          Sheet["D"*string(ind)] = y
#         end
#     end

# ind = 101
# xs = 0.0:4.0/(ind-1):4.0
# ys = 0.0:4.0/(ind-1):4.0
# zs = zeros(ind,ind)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         zs[i,j] = 𝑢(x,y)
#     end
# end

# fig = Figure()
# ax = Axis3(fig[1,1])
# surface!(ax,xs,ys,zs)
# fig

# xs = 0.0:0.4:4.0
# ys = 0.0:0.4:4.0
# zs = hcat([d[1],d[40:-1:32]...,0.0],
#           [d[5],d[41:49]...,0.0],
#           [d[6],d[50:58]...,d[30]],
#           [d[7],d[59:67]...,d[29]],
#           [d[8],d[68:76]...,d[28]],
#           [d[9],d[77:85]...,d[27]],
#           [d[10],d[86:94]...,d[26]],
#           [d[11],d[95:103]...,d[25]],
#           [d[12],d[104:112]...,d[24]],
#           [d[13],d[113:121]...,d[23]],
#           [d[2],d[14:22]...,d[3]])
# xs = zeros(nₚ)
# ys = zeros(nₚ)
# zs = zeros(nₚ)
# for (i,node) in enumerate(nodes)
#     xs[i] = node.x
#     ys[i] = node.y
#     zs[i] = node.d
# end

# fig


    