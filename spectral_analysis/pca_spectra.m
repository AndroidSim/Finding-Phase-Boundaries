function [RE,IE,IND,chisq] = pca_spectra(spectra,mean_center_flag,varargin)
% spectra = a (b x s) matrix, b = # of magnetic field point, s = # of samples
% spectra must be aligned and only contain the intensity values

[nb,ns] = size(spectra); % np = number of magnetic field points of each spectrum, ns = number of composition samples

if any(size(mean_center_flag) ~= 1)
    error('second argument must be a scalar');
else
    mean_center = logical(mean_center_flag);
end

if nargin > 2
    % smooth spectra along composition axis
    naa = length(varargin); % either 1 or 2
    
    if ~ischar(varargin{1})
        error('first argument of smoothing options must be a string specifying the smoothing method');
    else
        if any(strcmp(varargin{1},{'moving';'lowess'}))
            % 'moving' = moving average smoothing 
            % 'lowess' = local regression smoothing
            method = varargin{1};
        else
            error('first argument of smoothing options must be a string specifying the smoothing method');
        end
    end
    
    if naa > 1
        if ~isnumeric(varargin{2})
            error('second argument of smoothing options must be a number specifying the smoothing window size');
        else
            windowsize = varargin{2};
        end
    else
        switch method
            case 'moving'
                windowsize = 3; % default
            case 'lowess'
                windowsize = 5;
            otherwise
                windowsize = 3;
        end
    end
    
    for i = 1:nb
        spectra(i,:) = smooth(c,spectra(i,:),windowsize,method)';
    end
end

if mean_center
    spectra = spectra - repmat(mean(spectra,2),1,ns);
end

% svd to obtain principal components, eigenspectra, and eigenvalues
% [u,latent,pc] = svd(x,0);

[U,W,V] = svd(spectra,0);
evs = diag(W);
% U = eigenspectra
% W = eigenvalues
% V = eigenvectors

% determine number of principal components (n) necessary to describe data
for n = 1:ns-1
    RE(n) = sqrt(sum(evs(n+1:end))/(nb*(ns-n)));
    IE(n) = RE(n)*(sqrt(n/ns));
    IND(n) = RE(n)/((ns-n)^2);
    
    spectra_n = U(:,1:n)*W(1:n,1:n)*V(:,1:n)';
    dsq = (spectra - spectra_n).^2;
    % dsq = dsq./((0.5)^2);
    chisq_comp(n) = sum(sum(dsq));
    chisq_expd(n) = (nb-n)*(ns-n);
end

chisq = [chisq_comp;chisq_expd];

return    