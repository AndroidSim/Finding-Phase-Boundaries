function spectra = scale_spectra(spectra,scale)

[nb,ncol] = size(spectra);

if ncol < 2 || rem(ncol,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = ncol/2;

if ~isscalar(scale)
    error('scale must be a number');
end

% from asc2dat.m:
% max_mag=400; % has been used for years. I don't know why?
% spc(:,2)=spc(:,2)/max(spc(:,2))*max_mag;

% scale so that maximum intensity is the same for all spectra 
[maxI,imaxI] = max(spectra(:,2:2:end));
for s = 1:ns
    spectra(:,s*2) = (spectra(:,s*2)/maxI(s))*scale;
%     spectra(:,s*2) = (spectra(:,s*2)/spectra(imaxI(s),s*2))*scale;
end

return