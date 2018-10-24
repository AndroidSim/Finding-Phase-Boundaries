function a_spectra = deriv2abs(d_spectra)
% deriv2abs converts an derivative spectrum to its corresponding absorbance
% spectrum

[np,nc] = size(d_spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

for i = 1:ns
    a_spectra(:,i*2-1) = d_spectra(:,i*2-1);
    
    a_spectra(:,i*2) = cumsum(d_spectra(:,i*2));    
%     baseline correction
    d1 = d_spectra(1,i*2-1);
    a1 = a_spectra(1,i*2);
    d2 = d_spectra(end,i*2-1);
    a2 = a_spectra(end,i*2);
    A = (a2-a1)/(d2-d1); 
    B = a2-d2*A;
    a_spectra(:,i*2) = a_spectra(:,i*2)-(A*d_spectra(:,i*2-1)+B);
    
%     dH = zeros(size(d_spectra(:,i*2-1)));
%     dH(2:end) = diff(d_spectra(:,i*2-1));
%     a_spectra(:,i*2) = cumsum(prod([dH d_spectra(:,i*2)],2));
%     a_spectra(a_spectra(:,i*2) < 0,i*2) = 0;     
end

return