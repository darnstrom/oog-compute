## Init
include("oognorm.jl")

using ControlSystemsBase
using LinearAlgebra 

## Run experiments
# Quadruple tank (model taken from Teixeira 2021)
Ap = [-0.1068 0 0.0275 0;
      0 -0.0903 0 0.0258;
      0 0 -0.0275 0;
      0 0 0 -0.0258]
Bp = [0.0802 0;
      0 0.0807;
      0 0.1345;
      0.1337 0;]
Cm = [0.2 0 0 0; 0 0.2 0 0]
Q = diagm([40,40,20,20.0])
R = diagm([10,10])

Cp = [sqrt.(Q);zeros(2,4)]
Dp = [zeros(4,2);sqrt.(R)]

L = -lqr(Continuous,Ap,Bp,Q,R);
K = [0.0345 0.0150;
     0.0150 0.0414;
     0.0456 0.0633;
     0.0561 0.0515]

nx = size(Ap,1)
nu = size(Bp,2)  
ny = size(Cm,1)

E = [0,1e-9,1e-6,1e-3]
Ω=10 .^range(start=-4,stop=4,length=1000)
Ss =  [("unbounded", [2],[],[]),("ex7", [],[],[1])]

γs = zeros(length(Ω),length(E),length(Ss))
for (j,S) in enumerate(Ss)
    name,u_ind,y_ind,ext_ind = S 
    A = [Ap + Bp*L -Bp*L; 
         zeros(size(Ap)) Ap - K*Cm]

    Γu = I(nu)[:,u_ind]
    Γy = I(ny)[:,y_ind]
    Ep = I(nx)[:,ext_ind]

    B = [Bp*Γu zeros(size(Ap,1),size(Γy,2)) Ep;
         Bp*Γu K*Γy Ep]
    Cc = [Cp+Dp*L -Dp*L]  
    Dc = [Dp*Γu zeros(size(Cc,1),size(Γy,2)+size(Ep,2))]

    Cr = [zeros(size(Cm,1),size(Cp,2)) Cm]
    Dr = [zeros(size(Cr,1),size(Γu,2)) Γy zeros(size(Cr,1),size(Ep,2))]

    Gc, Gr  = tf(ss(A,B,Cc,Dc)), tf(ss(A,B,Cr,Dr))
    Gc0,Gr0 = dcgain(Gc),dcgain(Gr)
    println(sqrt(Gc0'*Gc0))
    println(sqrt(Gr0'*Gr0))
    for (k,ϵr) in enumerate(E)
        γs[:,k,j] = [oog(ω,A,B,Cc,Dc,Cr,Dr;ϵr) for ω in Ω]
    end

    # Write to file
    using DelimitedFiles
    open("result/reg_"*name*".dat"; write=true) do f
        write(f, "w e0 em9 em6 em3\n")
        writedlm(f,[Ω γs[:,:,j]])
    end
end

## Plotting
using Plots
pl1 = plot()
pl2 = plot()
for i in 1:4
    plot!(pl1,Ω,γs[:,i,2],xaxis=:log,ylims=(30,40))
    plot!(pl2,Ω,γs[:,i,1],xaxis=:log,yaxis=:log)
end
display(pl1)
display(pl2)
