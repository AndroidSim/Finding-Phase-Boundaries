function processed_spectra = spectra_preprocessing(spectra,varargin)
% spectra_preprocessing processes spectra to be used for the analysis of
% the variances of spectra due to the compositional differences of the
% bilayer.

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
center = varargin{1};
processed_spectra = spectra;

switch center
    case 'mean'
        % center the spectra about their mean
        for i = 1:np
            processed_spectra(i,2:2:end) = spectra(i,2:2:end)-mean(spectra(i,2:2:end));
        end
    case 'reference'
        ref_spectrum = varargin{2};
        [rnp,rnc] = size(ref_spectrum);

        if rnc ~= 2
            error('reference spectrum must contain 2 columns: [B-field absorbance_values]');
        end

        if rnp ~= np
            error('the number of points of the reference spectrum must equal the number of points of the processed spectra');
        end
        
        rns = rnc/2;
        
        % mean center the spectra with respect to a reference spectrum
        for i = 1:np
            processed_spectra(i,2:2:end) = spectra(i,2:2:end)-ref_spectrum(i,2));
        end
    otherwise
        error('no type of centering preprocessing was chosen');
end

return