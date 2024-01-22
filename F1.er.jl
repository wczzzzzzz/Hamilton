function (op::Operator{:∫q̇mpqkpdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    m = op.m
    kᶜ = op.kᶜ
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (B[i]*m*B[j] - N[i]*kᶜ*N[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫q̇mpqkpdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    ρ = op.ρ
    A = op.A
    EA = op.EA
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (Bₜ[i]*ρ*A*Bₜ[j] - Bₓ[i]*EA*Bₓ[j])*𝑤
            end
            f[I] += -(ρ*A*Bₜ[i])*N[i]*𝑤
        end
    end
end
# f[I] += -(ρ*A*Bₜ[i])*N[j] + ((N[i]*b) + N[j]*P)*𝑤
