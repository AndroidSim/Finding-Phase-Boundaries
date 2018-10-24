function S = get_spectrum(spectra,sn)
% get_spectrum:  retrieves a single spectrum from a set of spectra given by
% the spectrum number sn

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

if ns == 1
    S = spectra;
    return
end

ci = [sn*2-1 sn*2]
S = spectra(:,ci)

return