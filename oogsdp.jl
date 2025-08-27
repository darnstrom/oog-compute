## SDP 
using JuMP
using SparseArrays
using Clarabel

function oogsdp(A,B,Cr,Dr,Cc,Dc;ϵ=0, model=Model(Clarabel.Optimizer), isPSD=true, isDiagonal=false)
    n,m = size(B)

    if isDiagonal
        @variable(model,x[1:n])
        X = diagm(x)
    else
        if isPSD
            @variable(model, X[1:n, 1:n], PSD)
        else
            @variable(model, X[1:n, 1:n])
        end
    end
    @variable(model, γ≥0)

    G = [A'*X + X*A X*B; B'*X -ϵ*γ*I(m)]-γ*[Cr Dr]'*[Cr Dr] + [Cc Dc]'*[Cc Dc];

    @constraint(model, -sparse(G) in JuMP.PSDCone())
    @constraint(model, γ ≤ 1e5)
    @objective(model, Min, γ)

    optimize!(model)
    solution_summary(model)
    return sqrt(value(γ))
end
