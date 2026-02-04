function sigmabar(ω,A,B,C,D;a=0)
    if(ω ≠ Inf)
        M = ((im*ω+a)*I-A)\B
        T = C*M+D
    else
        T = D
    end
    λ = size(B,2) == 1 ? T : eigvals!(T'T);
    return real(sqrt(λ[end]))
end
function oog(P,ω;ϵr=0, a=0)
    return oog(ω,P.A,P.B1,P.C1,P.D11,P.Cr,P.Dr1;ϵr, a)
end
function oog(P,Ω::Matrix{<:Complex};ϵr=0, a=0)
    return [oog(P,ω;ϵr,a) for ω in Ω] 
end
function oog(P,Ω::Vector{<:Real};ϵr=0, a=0)
    return [oog(P,ω;ϵr,a) for ω in Ω] 
end

function oog(ω,A,B,Cp,Dp,Cr,Dr;ϵr=0, a=0, type=:naive)
    if(ω ≠ Inf)
        M = ((im*ω+a)*I-A)\B
        Tp, Tr = Cp*M+Dp, Cr*M+Dr
    else
        Tp, Tr = Dp, Dr
    end
    if(type==:naive)
        λ = size(B,2) == 1 ? Tp'*Tp/(Tr'*Tr+ϵr*I) : eigvals!(Tp'*Tp,Tr'*Tr+ϵr*I);
        return real(sqrt(λ[end]))
    else
        return maximum(svdvals!(Tp,[Tr;sqrt(ϵr)I]))
    end
end

function oognorm(P;ϵr=1e-12,tolimag=1e-5,ϵ=1e-5,a=0, max_iter=20)
    return oognorm(P.A,P.B1,P.C1,P.D11,P.Cr,P.Dr1;ϵr,tolimag,ϵ,a, max_iter)
end

#  ==== Algorithm 1 in paper  ====
function oognorm(A,B,Cp,Dp,Cr,Dr;ϵr=1e-12,tolimag=1e-5,ϵ=1e-5,a=0, max_iter=20, get_secondary=false)
    A = A-a*I;

    # Try ω = 0
    γl = oog(0.0,A,B,Cp,Dp,Cr,Dr;ϵr,a)/(1+2*ϵ)
    Ω_primary, Ω2= [0.0],[]; 
    if(isnan(γl))
        γl = -Inf
    end

    # Try ω = ∞
    sλ = oog(Inf,A,B,Cp,Dp,Cr,Dr;ϵr,a)
    if(sλ > γl*(1+2*ϵ))
        γl,Ω_primary = sλ,[Inf]
    end

    for i in 1:max_iter
        println("iter: $i | $(sqrt(γl))")
        γ = (1+2*ϵ)*γl

        H = hamiltonian(A,B,Cp,Dp,Cr,Dr,γ,ϵr)
        Ω = get_frequencies(H;tolimag)

        if(!isempty(Ω)) # Empty Ω => γopt ≥ γ
            γl = γ
        else
            if(get_secondary)
                Ht = hamiltonian(A,B,Cp,Dp,Cr,Dr,γ*0.95,ϵr)
                Ω2 = get_frequencies(Ht;tolimag)
            end
            return (1+ϵ)*γl,Ω_primary,Ω2
        end

        if(length(Ω) == 1)
            Ωm = [0.5*Ω[1], Ω[1], 1.5*Ω[1]]
        else
            Ωm = [0.5*(Ω[i+1] + Ω[i]) for i in 1:length(Ω)-1]
        end
        Ω_primary = Float64[] 
        for ω in Ωm 
            sλ = oog(ω,A,B,Cp,Dp,Cr,Dr;ϵr, a)
            γl = max(γl,sλ)
            if(sλ ≥ 0.995*γ)
                push!(Ω_primary,ω)
            end
        end
    end
    return (1+ϵ)*γl,Ω_primary,Ω2
end

function hamiltonian(A,B,C1,D1,Cr,Dr,γ,ϵ)
    iDD = inv(γ^2*(Dr'*Dr+ϵ*I)-D1'*D1)
    K = iDD*(D1'*C1-γ^2*Dr'*Cr)
    At = A+B*K
    H = [At -B*iDD*B'; 
         -γ^2*Cr'*(Cr+Dr*K)+C1'*(C1+D1*K) -At']
    return H
end

function get_frequencies(H; tolimag=1.5*1e-7)
        λs = eigvals(H) 
        Ω = [imag(λ) for λ in λs 
             if imag(λ) ≥ 0 && abs(real(λ))< tolimag*(1+abs(λ))]
        return sort!(Ω) 
end
