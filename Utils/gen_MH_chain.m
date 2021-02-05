function [MHchains]=gen_MH_chain(chain_size,burnin,num_restarts,v_list,N,evals,varargin)
    
    eff_size=chain_size/num_restarts;
    MHchains=zeros([max(size(v_list))*chain_size,N]); sigz=zeros([max(size(v_list)),1]);
    temp_chain=zeros([num_restarts*eff_size,N]);

    for ii=1:max(size(v_list))    

        for ll=1:num_restarts
            options=varargin;
            %options={'mps_psi', all_psi{ii},'angle',rot};
            burnchain=MetropolisHastings(burnin,evals,options{:});

            options = [{'s0',burnchain(end,:)} options];
            temp_chain((ll-1)*eff_size+1:(ll*eff_size),:)=...
                MetropolisHastings(eff_size,evals,options{:});
        end

        MHchains((ii-1)*chain_size+1:(ii*chain_size),:)=temp_chain;
        %sigz(ii)=mean(mean(MHchains((ii-1)*chain_size+1:(ii*chain_size),:)));
    end


end