% function Z=getZ(N,S)
%     D=2*S+1;
%     Z=zeros(D^N,N); 
%     ia=(0:D^N-1)';
%     for m=1:N 
%         ja=floor(ia/D);
%         Z(:,N+1-m)=S-ia+D*ja;
%         ia=ja;
%     end
% end

function Z=getZ(N,S)
% Modified code originally by my Advisor: Zhexuan Gong
%The linear-spin basis conversion. Tested to be the Fastest algorithm on matlab
%Z(:,m) is the diagonal elements of spin operator Z_m

D=2*S+1;
Z=zeros(D^N,N); 
ia=(0:D^N-1)';
for m=1:N 
    ja=floor(ia/D);
    Z(:,N+1-m)=S-ia+D*ja;
    ia=ja;
end

end
