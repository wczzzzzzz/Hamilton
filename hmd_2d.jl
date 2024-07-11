using  ApproxOperator, JuMP, Ipopt, CairoMakie, XLSX, LinearAlgebra

using GLMakie

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

model = Model(Ipopt.Optimizer)

include("import_hmd_test.jl")

ndiv= 10
# elements,nodes = import_hmd_Tri3("./msh/Non-uniform_"*string(ndiv)*".msh")
elements,nodes = import_hmd_Tri3("./msh/square_"*string(ndiv)*".msh")
# elements,nodes = import_hmd_Tri3("./msh/bar_"*string(ndiv)*".msh")
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
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
U(x) = α/π*(sin(π*(4-x/α)))
prescribe!(elements["Γ₁"],:𝑃=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
# prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->U(x))
prescribe!(elements["Γ₃"],:𝑃=>(x,y,z)->0.0)
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
       # Operator{:∫vgdΓ}(:α=>1.0),
    #    Operator{:L₂}(),
]

ops[1](elements["Ω"],k)
ops[2](elements["Γ₁"],f)
ops[2](elements["Γ₃"],f)
ops[3](elements["Γ₄"],f)
# ops[4](elements["Γ₁"],kᵅ,fᵅ)
ops[4](elements["Γ₂"],kᵅ,fᵅ)
ops[4](elements["Γ₃"],kᵝ,fᵝ)

# dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# d = dt[1:nₚ]
# δd = dt[nₚ+1:end]

# push!(nodes,:d=>d,:δd=>δd)

d = (k+kᵅ)\(f+fᵅ)
# d = k\f
push!(nodes,:d=>d)

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
# L₂ = ops[5](elements["Ωᵍ"])

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     Δ = d[i] - 𝑢(x,y)
#         index = [10,20,40,80]
#         XLSX.openxlsx("./excel/Non-uniform.xlsx", mode="rw") do xf
#         Sheet = xf[4]
#         ind = findfirst(n->n==ndiv,index)+i
#         Sheet["A"*string(ind)] = x
#         Sheet["B"*string(ind)] = y
#         Sheet["C"*string(ind)] = d₁
#         Sheet["D"*string(ind)] = Δ
#         # Sheet["E"*string(ind)] = log10(L₂)
#         # Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end

# push!(nodes,:d=>d₁)
# for (i,a) in enumerate(elements["Ω"])
#     node1 = a.𝓒[1]
#     node2 = a.𝓒[2]
#     node3 = a.𝓒[3]
#     x1 = node1.x
#     x2 = node2.x
#     x3 = node3.x
#     y1 = node1.y
#     y2 = node2.y
#     y3 = node3.y
#     d1 = node1.d
#     d2 = node2.d
#     d3 = node3.d
#     Δ1 = d1 - 𝑢(x1,y1)
#     Δ2 = d2 - 𝑢(x2,y2)
#     Δ3 = d3 - 𝑢(x3,y3)
#     XLSX.openxlsx("./excel/hmd_2d_n=10.xlsx", mode="rw") do xf
#         Sheet = xf[1]
#         ind = findfirst(n->n==ndiv,11)+6*(i-1)
#         Sheet["A"*string(ind+1)] = x1
#         Sheet["B"*string(ind+1)] = y1
#         Sheet["C"*string(ind+1)] = d1
#         Sheet["D"*string(ind+1)] = Δ1

#         Sheet["A"*string(ind+2)] = x2
#         Sheet["B"*string(ind+2)] = y2
#         Sheet["C"*string(ind+2)] = d2
#         Sheet["D"*string(ind+2)] = Δ2

#         Sheet["A"*string(ind+3)] = x3
#         Sheet["B"*string(ind+3)] = y3
#         Sheet["C"*string(ind+3)] = d3
#         Sheet["D"*string(ind+3)] = Δ3

#         Sheet["A"*string(ind+4)] = 0.5*(x1+x2)
#         Sheet["B"*string(ind+4)] = 0.5*(y1+y2)
#         Sheet["C"*string(ind+4)] = 0.5*(d1+d2)
#         Sheet["D"*string(ind+4)] = 0.5*(Δ1+Δ2)

#         Sheet["A"*string(ind+5)] = 0.5*(x2+x3)
#         Sheet["B"*string(ind+5)] = 0.5*(y2+y3)
#         Sheet["C"*string(ind+5)] = 0.5*(d2+d3)
#         Sheet["D"*string(ind+5)] = 0.5*(Δ2+Δ3)

#         Sheet["A"*string(ind+6)] = 0.5*(x3+x1)
#         Sheet["B"*string(ind+6)] = 0.5*(y3+y1)
#         Sheet["C"*string(ind+6)] = 0.5*(d3+d1)
#         Sheet["D"*string(ind+6)] = 0.5*(Δ3+Δ1)
#     end
# end

# push!(nodes,:d=>d)
# fig = Figure()
# Axis(fig[1, 1])
# xs = [node.x for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# ys = [node.d for node in nodes[[36,45,54,63,72,81,90,99,108,117,18]]]
# lines!(xs,ys, color = :blue)

# fig

# ind = 101
# xs = 0.0:4.0/(ind-1):4.0
# ys = 0.0:4.0/(ind-1):4.0
# zs = zeros(ind,ind)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         zs[i,j] = 𝑢(x,y)
#     end
# end

fig = Figure()
ax = Axis3(fig[1,1])
# fig

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
    # δds[i] = node.δd
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
meshscatter!(ax,xs,ys,ds,color=ds,markersize = 0.1)
# meshscatter!(ax,xs,ys,δds,color=δds,markersize = 0.1)
fig


    