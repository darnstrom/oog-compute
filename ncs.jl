using Graphs
using LinearAlgebra
using Mosek, MosekTools
using PProf, Profile
include("oognorm.jl")
include("oogsdp.jl")
## Experiments
import Random
Random.seed!(1234)
Ns = 50:50:1000
nrun = 50
ts = fill(NaN,3,length(Ns),nrun)
γs = fill(NaN,3,length(Ns),nrun)
for (i,N) in enumerate(Ns)
    for j in 1:nrun
        nd,na = Int(N/50),Int(N/50)
        g = SimpleDiGraph(N, Int(round(1.0*N)))
        L = laplacian_matrix(g)
        cc = connected_components(g)
        # Connected graph
        for i = 1:length(cc)-1
            add_edge!(g,rand(cc[i]),rand(cc[i+1]))
        end

        # Dynamics 
        θ = 1.0
        A = float.(adjacency_matrix(g))
        for e in edges(g)
            A[src(e),dst(e)] = 1+2e-1*rand()
        end
        A = Graphs.LinAlg.symmetrize(A)
        A -= diagm(θ .+ sum(A,dims=2)[:]) # θ for self loop

        # Performance metric
        Cp = diagm(ones(N))
        Dp = zeros(N,na) 

        # Detector
        v = collect(1:N)
        d_ids = [popat!(v,rand(1:length(v))) for i in 1:nd]
        Cr = I(N)[d_ids,:]
        Dr = zeros(nd,na)

        # Attack signal
        v = collect(1:N)
        a_ids = [popat!(v,rand(1:length(v))) for i in 1:na]
        B = I(N)[:,a_ids]

        # Compute
        ts[1,i,j] = @elapsed γs[1,i,j],ω,_ = oognorm(A,B,Cp,Dp,Cr,Dr;ϵr=1e-5)

        model = Model(Clarabel.Optimizer)
        γs[2,i,j] = oogsdp(A,B,Cr,Dr,Cp,Dp;model,ϵ=1e-5,isPSD=false,isDiagonal=true)
        ts[2,i,j] = solve_time(model)
        if(N <= 100)
            model = Model(Mosek.Optimizer)
            γs[3,i,j] = oogsdp(A,B,Cr,Dr,Cp,Dp;model,ϵ=1e-5,isPSD=false,isDiagonal=false)
            ts[3,i,j] = solve_time(model)
        end
    end
end
## Collect

using Statistics
# ts 
tmax = maximum(ts,dims=3)[:,:,1]
tmin = minimum(ts,dims=3)[:,:,1]
tavg = mean(ts,dims=3)[:,:,1]
tmed = median(ts,dims=3)[:,:,1]
# correct
γham = γs[1,:,:][:]
γclara = γs[2,:,:][:]
diffs = abs.(γclara-γham)./γham

## Writing to file
using DelimitedFiles
for (id,name) in [(1,"hamiltonian"),(2,"clarabel"),(3,"mosek")]
    open("result/nwc_"*name*".dat"; write=true) do f
        write(f, "nx tmin tavg tmed tmax\n")
        writedlm(f,[Ns tmin[id,:] tavg[id,:] tmed[id,:] tmax[id,:]])
    end
end
open("result/nwc_diff.dat"; write=true) do f
    write(f, "errors\n")
    writedlm(f,[diffs;;])
end
## Try graphs which are connected more
N = 500
nd,na = Int(N/50),Int(N/50)
g = SimpleDiGraph(N, Int(round(1.5*N)))
L = laplacian_matrix(g)
cc = connected_components(g)
# Connected graph
for i = 1:length(cc)-1
    add_edge!(g,rand(cc[i]),rand(cc[i+1]))
end

# Dynamics 
θ = 1.0
A = float.(adjacency_matrix(g))
for e in edges(g)
    A[src(e),dst(e)] = 1+2e-1*rand()
end
A = Graphs.LinAlg.symmetrize(A)
A -= diagm(θ .+ sum(A,dims=2)[:]) # θ for self loop

# Performance metric
Cp = diagm(ones(N))
Dp = zeros(N,na) 

# Detector
v = collect(1:N)
d_ids = [popat!(v,rand(1:length(v))) for i in 1:nd]
Cr = I(N)[d_ids,:]
Dr = zeros(nd,na)

# Attack signal
v = collect(1:N)
a_ids = [popat!(v,rand(1:length(v))) for i in 1:na]
B = I(N)[:,a_ids]

model = Model(Clarabel.Optimizer)
γc = oogsdp(A,B,Cr,Dr,Cp,Dp;model,ϵ=1e-5,isPSD=false,isDiagonal=true)
tc = solve_time(model)
