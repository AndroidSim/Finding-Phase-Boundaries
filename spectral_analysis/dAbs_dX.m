function F = dAbs_dX(spectra,X,v)
% dAbs_dX: calculates the delta absorbance over delta composition of a set
% of spectra and their compositions
%
% *assumption*: spectra are aligned 
%
% spectra = (number points x number of spectra) matrix, 
%   spectra just contains absorbances values
% X = (v x number of spectra) vector
% v = the number of composition variables

[np,ns] = size(spectra);

%if ns ~= 2
%    error('number of spectra must = 2');
%end

[m,n] = size(X);

if m ~= v || n ~= ns
    error('X must be a row vector of size number variables x number of spectra');
end

F = zeros(np,ns-1);
for i = 1:np
    F(i,:) = diff(spectra(i,:),1,2)./diff(X); % == (spectra(:,2)-spectra(:,1))/(X(2)-X(1));
% [FX,FY] = gradient(spectra,0.05,diff(X)); F = FX;
end