function norm_data(data_file)
% norm_data normalizes the spectrum data given in the file data_file and
% saves the normalized spectrum to a file

if (exist(data_file,'file') == 0)
    error('data file does not exist');
end

file = which(data_file);
[pathstr,filename,ext,versn] = fileparts(file);

if ~strcmp(ext,'.dat')
    error('data file must be a .dat file');
end

spc = load([filename,ext]);
[nr,nc] = size(spc);
nf = nc/2;

if nf > 1
    norm_dspc = zeros(size(spc));
    norm_aspc = zeros(size(spc));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% the following code was copied and slightly modified from Yun-Wei's
    %% program DPSpectra_ct
    
    % creating normalized difference (normal) and absorbance spectrum
    sm=1;sn=1024;b0=3326;
    
    for j=2:2:nc
        spca1=spc;spca1(1,j)=0;
        
        for i=sm+1:sn
            spca1(i,j)=spca1(i-1,j)+spc(i,j);%ABS spectrum
        end
        
        %baseline correction
        a1=spc(1,j-1);b1=spca1(1,j);a2=spc(sn,j-1);b2=spca1(sn,j);
        A=(b2-b1)/(a2-a1); B=b2-a2*A;
        spca1(sm:sn,j)=spca1(sm:sn,j)-(A*spc(sm:sn,j-1)+B);

        spca1(sm:sn,j)=spca1(sm:sn,j)/sum(spca1(sm:sn,j))/mean...
        (spca1((sm+1):sn,j-1)-spca1(1:(sn-1),j-1));%Normalize ABS spectrum
        
        % just for difference spectra
        for i=sm+1:sn
            spc(i,j)=spca1(i,j)-spca1(i-1,j);
        end
        spc(1,j)=0;
        [c1,i1]=max(spc(:,j));
        displac=spc(i1,j-1)-b0;
        spc(:,j-1)=spc(:,j-1)-displac;
        spc1n(:,1)=spc(:,j-1);
        spc1n(:,2)=spc(:,j);
        norm_dspc(:,j-1:j) = spc1n(:,1:2);
    
        % just for absorbance spectra
        [c1,i1]=max(spca1(:,j));
        displac=spca1(i1,j-1)-b0;
        spca1(:,j-1)=spca1(:,j-1)-displac;
        spc1n(:,2)=spca1(:,j);
        spc1n(:,1)=spca1(:,j-1);
        norm_aspc(:,j-1:j) = spc1n(:,1:2);
    end 
    
    clear spcln;
else
    % creating normalized difference (normal) spectrum 
    norm_dspc = DPSpectra_ct(filename);  

    % creating normalized absorbance spectrum
    norm_aspc = DPSpectra_ct(filename,1);
end

clear spc;
sd_filename = ['normd_',filename,ext];
sa_filename = ['norma_',filename,ext];
format = [repmat('%8.3f %10.5f ',1,nf-1), '%8.3f %10.5f\n'];
dfid = fopen(sd_filename,'w');
afid = fopen(sa_filename,'w');
fprintf(dfid,format,norm_dspc');
fclose(dfid);
fprintf(afid,format,norm_aspc');
fclose(afid);