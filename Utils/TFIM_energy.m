function exact_energy=TFIM_energy(N,v_list)

exact_energy=zeros([max(size(v_list)),1]);
for i=1:max(size(v_list))
    h=v_list(i);
    if mod(N,2)==0
        qlist=(2*pi/N)*(-0.5*(N-1):(0.5*N));
    else
        qlist=(2*pi/N)*((-0.5*N):(0.5*(N+1)));
    end
    exact=0;
    for x=1:(N-1)
        wq=sqrt(1+2*h*cos(qlist(x))+h^2);
        exact=exact-wq;
    end
    exact_energy(i)=exact;
end

end