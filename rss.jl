using LinearAlgebra
# Generate a random statble state-space model
function rss(nx,ny,nu; ndims = nothing)
    # Generate A
    A = zeros(nx,nx);

    if(isnothing(ndims) || length(ndims)!=4)
        nint = (rand()<0.10)+sum(rand(nx-1).<0.01);
        nint = 0; # No integrals...
        nrepeated = Int(floor(sum(rand(nx-nint).<0.05)/2));
        ncomplex = Int(floor(sum(rand(nx-nint-2*nrepeated).<0.5)/2));
        nreal = nx-nint-2*nrepeated-2*ncomplex;
    else
        nint,nrepeated,ncomplex,nreal = ndims
        if nint+2*nrepeated+2*ncomplex+nreal != nx  
            @error("inconsistent specification of ndims")
        end
    end

    # Poles
    repeated_poles = -exp.(randn(nrepeated));
    complex_re = -exp.(randn(ncomplex));
    complex_im = 3 * exp.(randn(ncomplex));

    # Copmlex poles
    for i in 1:ncomplex
        A[2*i-1:2*i,2*i-1:2*i] = [ complex_re[i] complex_im[i];
                                  -complex_im[i] complex_re[i]];
    end

    # Integrator + repeated poles + real poles
    A[2*ncomplex+1:nx,2*ncomplex+1:nx] = diagm([zeros(nint);
                                              repeated_poles;
                                              repeated_poles;
                                              -exp.(randn(nreal))]);
    # Transform 
    #T = qr(randn(nx,nx)).Q;
    T = randn(nx,nx)
    A = T\A*T;
    # Generate B,C,D
    B,C,D = randn(nx,nu), randn(ny,nx), randn(ny,nu)
    # TODO: Add some sparsity
    return A,B,C,D
end
