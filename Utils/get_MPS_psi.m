function [all_psi, v_list,energy_list] = get_MPS_psi(model,L)
% this function generates and reshapes the psi matrix from the OSMPS
% generated .log files. It returns the psi for each variable in the
% variable list v_list in the form of a cell. 

% the output directory of the MPS sims
mps_dir=['/home/alex/Documents/QML_Research/NetKet/'...
    'Output/MPS_Results/' model '/' model '_OUTPUTS' num2str(L) '/'];

% this function gets a list of all .log files in the directory
out_files=glob([mps_dir '**.log']);

all_psi={}; v_list=[]; energy_list=[];
for lol=1:max(size(out_files)) % loops through each file
    file=out_files{lol};
    g_s= strfind(file,'g_');
    v_list(end+1)=str2num(file(g_s(end)+2:strfind(file,'.log')-1));
    
    % read the file as a str
    dat=fileread(file);
    
    % find the energy
    energy_list(end+1)=str2num(dat(strfind(dat,'Energy: ')+8:strfind(dat,'with variance=')-1));

    % find the beggining and end of the list of numbers that is Psi
    start=strfind(dat,'Psi :'); old=false;
    if isempty(start) % for older dat file versions
        start=strfind(dat,'Psi:'); old=true;
    end
    ending = strfind(dat,'time taken');
    if contains(dat(start+5:ending-1),')') %means the output is complex
        psi_list=str2num(erase(erase(dat(start+5:ending-1),'('),')'));
        iscomplex=true;
    else
        psi_list=str2num(dat(start+5:ending-1)); iscomplex=false;
    end
    
    % Similar idea to extract the matrix sizes of the A matrices that make up
    % psi
    A_dim1= str2num(dat(strfind(dat,'A%dl(1):')+8:strfind(dat,'A%dl(2)')-1));
    A_dim2= str2num(dat(strfind(dat,'A%dl(2):')+8:strfind(dat,'A%dl(3)')-1));
    if old
        A_dim3= str2num(dat(strfind(dat,'A%dl(3):')+8:strfind(dat,'Psi:')-1));
    else
        A_dim3= str2num(dat(strfind(dat,'A%dl(3):')+8:strfind(dat,'Psi :')-1));
    % dim2 is always special in the output, it is the size of the local hilbert space
    end
    
    % reshaping 
    psi={}; count=1;
  
    for ii=1:L
        AL = zeros([A_dim1(ii),A_dim2(ii),A_dim3(ii)]);
        for j1=1:A_dim1(ii)
            for j2=1:A_dim2(ii)
                for j3=1:A_dim3(ii)
                if iscomplex
                    AL(j1,j2,j3)= complex(psi_list(count),psi_list(count+1));
                    count=count+2;
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
    
end



end 