function d_spectra = abs2deriv(a_spectra)
% abs2deriv converts an absorbance spectrum to its corresponding derivative
% spectrum

[np,nc] = size(a_spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
d_spectra = zeros(size(a_spectra));

for i = 1:ns
    % convert absorbance spectrum to derivative spectrum
    d_spectra(2:end,i*2) = diff(a_spectra(:,i*2));
    d_spectra(:,i*2-1) = a_spectra(:,i*2-1);
end

return