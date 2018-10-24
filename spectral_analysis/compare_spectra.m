function delta_spec = compare_spectra(spectra,all2one_flag)
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
    else
        all2one = false;
    end
end
 
if all2one
    % align spectra about central peak with respect to first spectrum of
    % spectra
    alined_spec = align_spectra(spectra);

    % calculate difference spectra (delta_spec)
    delta_spec = ones(size(spectra))*sqrt(NaN);
    delta_spec(:,1:2) = spectra(:,1:2);

    for i=2:ns
        ind = find(alined_spec(:,i*2));
        delta_spec(ind,i*2) = spectra(ind,i*2)-spectra(ind,2); % difference spectrum
        delta_spec(ind,i*2-1) = spectra(ind,i*2-1);
    end
else
    % calculate difference spectra (delta_spec)
    delta_spec = ones(size(spectra))*sqrt(NaN);

    for i=1:ns-1
        % align pair of spectrum
        alined_spec = align_spectra(spectra(:,i*2-1:i*2+2));
        % find common indices of aligned spectra
        ind1 = find(alined_spec(:,2));
        ind2 = find(alined_spec(:,4));
        ind = intersect(ind1,ind2);
        delta_spec(ind,i*2) = alined_spec(ind,4)-alined_spec(ind,2); % difference spectrum
        delta_spec(ind,i*2-1) = alined_spec(ind,3);
    end

    alined_spec = align_spectra([spectra(:,end-1:end) spectra(:,1:2)]);
    ind1 = find(alined_spec(:,2));
    ind2 = find(alined_spec(:,4));
    ind = intersect(ind1,ind2);
    delta_spec(ind,end) = alined_spec(ind,4)-alined_spec(ind,2); % difference spectrum
    delta_spec(ind,end-1) = alined_spec(ind,1);
end

% calculate the taxi-cab norm (norm of order 1) of difference spectra
% norm(spectra(:,i*2)-spectra(:,i*2+2),1) == sum(abs(spectra(:,i*2)-spectra(:,i*2+2)))
%for i=1:ns
%    norm_one(i) = norm(delta_spec(:,i*2),1); 
%end

%plot(norm_one,'.k','markersize',20)