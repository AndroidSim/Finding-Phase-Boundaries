function spc1n = DPSpectra_ct(filnm1,abs_flag);
%%Experimental spectrum Data Pre-processing program.
%%For generating C_spectra for the input of LLS_PT2_ct() program.
%%
%%==>   modified_spectrum = DPSpectra_ct(file_name_without_ext);
%%==>   modified_ABS_spectrum = DPSpectra_ct(file_name_without_ext,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%The input file must be a .dat file, which is a (1024) by 2 matrix.
%%
%%
%%On return: 
%%  1. A modified (Normanized, baseline-corrected and shifted) matrix (1024*2).
%%  2. The highest peak is aligned to 3326 G.
%%
%% by YWC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sm=1;sn=1024;b0=3326;
%sm2=420;sn2=579;
ext_dat='.dat';
filnm1=[filnm1,num2str(ext_dat)];
fid=fopen(filnm1);
spc = fscanf(fid,'%g',[2 inf]); % It has two rows now.
spc = spc';
fclose(fid);
spca1=spc;spca1(1,2)=0;
for i=sm+1:sn
    spca1(i,2)=spca1(i-1,2)+spc(i,2);%ABS spectrum
end
%baseline correction
a1=spc(1,1);b1=spca1(1,2);a2=spc(sn,1);b2=spca1(sn,2);
A=(b2-b1)/(a2-a1); B=b2-a2*A;
spca1(sm:sn,2)=spca1(sm:sn,2)-(A*spc(sm:sn,1)+B);
%
spca1(sm:sn,2)=spca1(sm:sn,2)/sum(spca1(sm:sn,2))/mean...
    (spca1((sm+1):sn,1)-spca1(1:(sn-1),1));%Normalize ABS spectrum
%
%[c1,i1]=max(spca1(:,2));
%displac=spca1(i1,1)-b0;
%spca1(:,1)=spca1(:,1)-displac;
if (nargin==1)
    for i=sm+1:sn
        spc(i,2)=spca1(i,2)-spca1(i-1,2);
    end
    spc(1,2)=0;
    [c1,i1]=max(spc(:,2));
    displac=spc(i1,1)-b0;
    spc(:,1)=spc(:,1)-displac;
    spc1n(:,1)=spc(:,1);spc1n(:,2)=spc(:,2);
    %
elseif (nargin==2 & abs_flag==1)
    [c1,i1]=max(spca1(:,2));
    displac=spca1(i1,1)-b0;
    spca1(:,1)=spca1(:,1)-displac;
    spc1n(:,2)=spca1(:,2);
    spc1n(:,1)=spca1(:,1);
end