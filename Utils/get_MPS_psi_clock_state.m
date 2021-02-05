function [dict] = get_MPS_psi_clock_state(model,L) % returns a dictionary
% this function generates and reshapes the psi matrix from the OSMPS
% generated .log files. It returns the psi for each variable in the
% variable list v_list in the form of a cell. It also returns other MPS 
% simulation details such as bond entanglement entropy, site entropy
% and whether the simulation converged and with what variance.

% the output directory of the MPS sims
mps_dir=['/home/alex/Documents/QML_Research/NetKet/Output/MPS_Results/'...
    model '/' model '_OUTPUTS' num2str(L)];

% this function gets a list of all .log files in the directory
out_files=glob([mps_dir '/**.log']); 
nfiles=max(size(out_files));

all_psi={}; v_list=[]; theta_list=[]; phi_list=[]; energy_list=[]; 
bondentropylist=zeros([nfiles,L+1]); siteentropylist=zeros([nfiles,L]);
avgtau=[]; avgsigma=[];converged=[];
for lol=1:nfiles % loops through each file
    tic
    file=out_files{lol};
    g_s= strfind(file,'g_');
    if contains(file,'theta') % for getting the input parameters from each file name
        theta_s=strfind(file,'theta_');  shift=max(size('g_'));
        v_list(end+1)=str2num(file((g_s(end)+shift):(theta_s(1)-1)));

        phi_s=strfind(file,'phi_'); shift=max(size('theta_'));
        theta_list(end+1)=str2num(file((theta_s(end)+shift):(phi_s(1)-1)));

        shift=max(size('phi_'));
        phi_list(end+1)=str2num(file(phi_s(end)+shift:strfind(file,'.log')-1));
    end
    
    % read the file as a str
    dat=fileread(file);

    % find the energy
    energy_list(end+1)=str2num(dat(strfind(dat,'Energy: ')...
        +8:strfind(dat,'with variance=')-1));
    
    converged(end+1)=boolean(dat(strfind(dat,'converged: ')+11)); %T if converged
    bondentropylist(lol,:)=str2num(dat(strfind(dat,'BondEntropy: ')+12:...
            strfind(dat,'avgtau')-1));
    avgtau(end+1)=str2num(dat(strfind(dat,'avgtau: ')+8:...
            strfind(dat,'avgsigma')-1));
    avgsigma(end+1)=str2num(dat(strfind(dat,'avgsigma: ')+10:...
            strfind(dat,'SiteEntropy')-1));
    siteentropylist(lol,:)=str2num(dat(strfind(dat,'SiteEntropy: ')+12:...
            strfind(dat,'converged:')-1));
    
    % find the beggining and end of the list of numbers that is Psi
    start=strfind(dat,'Psi :'); 

    ending = strfind(dat,'time taken');
    if contains(dat(start+5:ending-1),')') %means the output is complex
        psicell=strsplit(erase(dat(start+5:ending-1),'('),')');
        psi_list=[cellfun(@str2num,psicell(1:end-1),'un',0)];
        % psi_list=str2num(erase(erase(dat(start+5:ending-1),'('),')')); wasn't working
        iscomplex=true;
    else
        psi_list=str2num(dat(start+5:ending-1)); iscomplex=false;
    end
    
    % Similar idea to extract the matrix sizes of the A matrices that make up
    % psi
%     A_dim1= str2num(dat(strfind(dat,'A%dl(1):')+8:strfind(dat,'A%dl(2)')-1));
%     A_dim2= str2num(dat(strfind(dat,'A%dl(2):')+8:strfind(dat,'A%dl(3)')-1));
%     A_dim3= str2num(dat(strfind(dat,'A%dl(3):')+8:strfind(dat,'Psi :')-1));
% previously working code... (pre mio generated outputs)
    
    A1string=strsplit(dat(strfind(dat,'A%dl(1):')+8:strfind(dat,'A%dl(2)')-1));
    A_dim1=cell2mat([cellfun(@str2num,A1string(1:end),'un',0)]);
    
    % dim2 is always special in the output, it is the size of the local hilbert space
    A2string=strsplit(dat(strfind(dat,'A%dl(2):')+8:strfind(dat,'A%dl(3)')-1));
    A_dim2=cell2mat([cellfun(@str2num,A2string(1:end),'un',0)]);
    
    A3string=strsplit(dat(strfind(dat,'A%dl(3):')+8:strfind(dat,'Psi :')-1));
    A_dim3=cell2mat([cellfun(@str2num,A3string(1:end),'un',0)]);
    
    % reshaping 
    psi={}; count=1;
    for ii=1:L
        AL = zeros([A_dim1(ii),A_dim2(ii),A_dim3(ii)]);
        for j1=1:A_dim1(ii)
            for j2=1:A_dim2(ii)
                for j3=1:A_dim3(ii)
                if iscomplex
                    AL(j1,j2,j3)= complex(psi_list{count}(1),psi_list{count}(2));
                    count=count+1; % count=count+2; for previously working version
                else
                    AL(j1,j2,j3)= psi_list(count);
                    count=count+1;
                end
                end
            end
        end
        psi{end+1}=AL; 
    end
    
    all_psi{end+1}=psi;
    toc
end

% map all of the outputs to a dictionary for convient retrieval
dict=containers.Map;
dict('psi_list')=all_psi;
dict('v_list')=v_list;
dict('theta_list')=theta_list;
dict('phi_list')=phi_list;
dict('energy_list')=energy_list;
dict('bondentropy_list')=bondentropylist;
dict('siteentropy_list')=siteentropylist;
dict('avgtau_list')=avgtau;
dict('avgsigma_list')=avgsigma;
dict('converge_list')=converged;

end 