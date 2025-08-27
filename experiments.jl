## Init 
include("oognorm.jl")
include("oogsdp.jl")
include("rss.jl")

using ControlSystemsBase
using LaTeXStrings
using Mosek, MosekTools

## Run experiment
import Random
Random.seed!(1234)

ns = collect(1:10)
nruns = 100
ϵ = 1e-8
times = fill(NaN, 3,length(ns),nruns)
γs = zeros(3,length(ns),nruns)
rel_error= zeros(3,length(ns),nruns)
for n in ns
    for nrun in 1:nruns
        nx,ny,nu = 5*n,n,n
        A,B,Cc,Dc = rss(nx,ny,nu) 
        Cr,Dr = randn(ny,nx),randn(ny,nu)
        Gc, Gr  = tf(ss(A,B,Cc,Dc)), tf(ss(A,B,Cr,Dr))
        Dc[:] .= 0 ; Dr .= 0;
        P = (A=A,B1=B,C1=Cc,D11=Dc, Cr=Cr, Dr1 = Dr)
        tmy = @elapsed  γmy,om = oognorm(P;ϵr=ϵ,tolimag=1e-6)
        # ref
        Ω=10 .^range(start=-1,stop=4,length=10000) ∪ om
        tgrid = @elapsed γgrid = maximum(oog(P,Ω;ϵr=ϵ))

        times[1,n,nrun] = tmy
        γs[1,n,nrun] = γmy
        rel_error[1,n,nrun] = abs(γgrid-γmy)/γgrid

        times[2,n,nrun] = tgrid
        γs[2,n,nrun] = γgrid

        model = Model(Mosek.Optimizer)
        set_silent(model)
        γmosek = oogsdp(A,B,Cr,Dr,Cc,Dc;model,ϵ=ϵ,isPSD=false)
        tmosek= solve_time(model)
        times[3,n,nrun] = tmosek
        γs[3,n,nrun] = γmosek 
        rel_error[3,n,nrun] = abs(γgrid-γmosek)/γgrid
    end
end

## Collect

using Statistics
# times 
tmax = maximum(times,dims=3)[:,:,1]
tmin = minimum(times,dims=3)[:,:,1]
tavg = mean(times,dims=3)[:,:,1]
tmed = median(times,dims=3)[:,:,1]
# correct

acc = sum(rel_error .< 0.05,dims=3)[:,:,1]./nruns

## Plotting
plot(tavg[1,:])
plot!(tavg[3,:],yaxis=:log)

## Writing to file
using DelimitedFiles
for (id,name) in [(1,"hamiltonian"),(3,"mosek")]
    open("result/"*name*".dat"; write=true) do f
        write(f, "nx tmin tavg tmed tmax accuracy\n")
        writedlm(f,[5*ns tmin[id,:] tavg[id,:] tmed[id,:] tmax[id,:] acc[id,:]])
    end
end
