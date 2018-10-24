function n_spec = normal(spectra,aord)
% normal normalizes each spectrum in spectra
% so that the area under each spectrum = 1
%
% aord = 'a' for absorbance or 'd' for derivative

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
n_spec = zeros(size(spectra));

if aord == 'a'
    for i = 1:ns
        H_field = spectra(:,i*2-1);
        specP = spectra(:,i*2);
    
        % calculate area A
        A = spectral_prop([H_field specP],aord,'area');
        
        n_spec(:,i*2-1) = H_field;
        n_spec(:,i*2) = specP./A;
    end
end

if aord == 'd'
    for i = 1:ns
        H_field = spectra(:,i*2-1);
        specP = spectra(:,i*2);
    
        % calculate area A
        A = spectral_prop([H_field specP],aord,'area');
        
        n_spec(:,i*2-1) = H_field;
        n_spec(:,i*2) = specP./A;
    end
end

return