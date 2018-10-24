function spectra = difference_spectra(spectra,nofone)
% compare_spectra either compares 2 spectra or, if more than two spectra, 
% compare_spectra either compares all spectra to the first spectra if all2one_flag == 1
% or each spectrum is compared to the one before it.

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
all2one = false;

if ns > 2
    if nargin == 2
        all2one = true;
        
        if ~all(size(nofone) == 1)
            error('second argument must be a scalar');
        else
            if ~any(nofone == [1:ns])
                error('second argument must specify the spectrum number [1:number of spectra]');
            end 
        end
    end
end
 
% align spectra 
spectra = align_spectra(spectra);

% calculate difference spectra (delta_spec)
delta_spectra = spectra(:,1:end-2);

if all2one
    for s = 1:ns-1
        % ind = find(all(~isnan(alined_spec(:,i*2-1:i*2)),2));
        delta_spectra(:,s*2) = spectra(:,s*2)-spectra(:,nofone*2);
    end
else
    values = spectra(:,2:2:end);
    values = diff(values,1,2);
    
    for s = 1:ns-1
        % find common indices of aligned spectra
        % ind = find(all(~isnan(alined_spec),2));
        % ind1 = find(alined_spec(:,2));
        % ind2 = find(alined_spec(:,4));
        % ind = intersect(ind1,ind2);
        delta_spectra(:,s*2) = values(:,s); 
    end
end

spectra = delta_spectra;

return