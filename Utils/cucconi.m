% Greatly influenced by the R function written by Theo Pepler
function C=cucconi(x,y)

    m=numel(x); n=numel(y);
    N=m+n;
    dat=[x ; y];
    S=rank(dat(m+1:N));
    denom=sqrt(m*n*(N+1)*(2*N+1)*(8*N+11)/5);
    U=(6*sum(S.^2) - n*(N+1)*(2*N+1))/denom;
    V= (6*sum((N+1-S).^2) -n*(N+1)*(2*N+1))/denom;
    rho=(2*(N^2 -4))/((2*N+1)*(8*N+11))-1;
    
    C=(U^2+V^2-2*rho*U*V)/(2*(1-rho^2));
    
end