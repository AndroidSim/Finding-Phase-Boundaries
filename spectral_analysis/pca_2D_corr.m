function [varargout] = pca_2D_corr(spectra,c,reference,varargin)
% spectra = a (b x s) matrix, b = # of magnetic field point, s = # of samples
% spectra must be aligned and only contain the intensity values
% c = a (1 x s) vector of the compositions of the samples
% reference = the reference spectrum to be subtracted, if reference is not a vector,
% then no reference is subtracted
% varargin is for specifying smoothing, if needed.  see "smooth" matlab fxn

[nb,ns] = size(spectra); % np = number of magnetic field points of each spectrum, ns = number of composition samples

if all(size(c) == 1) | all(size(c) > 1) | ndims(c) > 2
    error('second argument must be a vector');
else
    nc = length(c); % nc = number of compositions, should equal ns
end

if ~isequal(nc,ns)
    error('the number of compositions must equal the number of samples');
end

if nargin > 3
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

% subtract reference spectrum (ie mean spectrum, single spectrum, etc.)
if any(size(reference) == 1) & ndims(reference) < 3
    if ~all(size(reference) == 1) % ie not scalar
        if ~isequal(nb,length(reference))
            error('number points of reference should equal number points of spectra');
        end
        
        r = size(reference);
        
        if r(1) == 1
            reference = reference';
            spectra = spectra - repmat(reference,1,ns);
        else
            spectra = spectra - repmat(reference,1,ns);
        end
    end
end

% calculate synchronous covariance matrices
cc_scov_mat = spectra' * spectra / (nb-1); % (c x b) * (b x c)
bb_scov_mat = spectra * spectra' / (nc-1); % (b x c) * (c x b)

% calculate asynchronous covariance matrices
% calculate H, the Hilbert-Noda matrix.  two forms: Hc for composition
% spectra and Hb for magnetic spectra
for m = 1:nb
    for n = 1:nb
        if m == n
            Hb(m,n) = 0;
        else
            Hb(m,n) = 1/(n-m);
        end
    end
end

for m = 1:nc
    for n = 1:nc
        if m == n
            Hc(m,n) = 0;
        else
            Hc(m,n) = 1/(n-m);
        end
    end
end

cc_acov_mat = spectra' * (Hb * spectra) / (nb-1); % (c x b) * [(b x b) * (b x c)]
bb_acov_mat = spectra * (Hc * spectra') / (nc-1); % (b x c) * [(c x c) * (c x b)]

% post filter for covariance matrices to set very small covariances to zero
bb_scov_max = max(max(bb_scov_mat));
bb_scov_min = min(min(bb_scov_mat));
bb_acov_max = max(max(bb_acov_mat));
bb_acov_min = min(min(bb_acov_mat));
cc_scov_max = max(max(cc_scov_mat));
cc_scov_min = min(min(cc_scov_mat));
cc_acov_max = max(max(cc_acov_mat));
cc_acov_min = min(min(cc_acov_mat));

for m = 1:nb
    for n = 1:nb
        if abs(bb_scov_mat(m,n)) < (bb_scov_max-bb_scov_min)/100
            bb_scov_mat(m,n) = 0;
        end
        if abs(bb_acov_mat(m,n)) < (bb_acov_max-bb_acov_min)/100
            bb_acov_mat(m,n) = 0;
        end
    end
end

for m = 1:nc
    for n = 1:nc
        if abs(cc_scov_mat(m,n)) < (cc_scov_max-cc_scov_min)/100
            cc_scov_mat(m,n) = 0;
        end
        if abs(cc_acov_mat(m,n)) < (cc_acov_max-cc_acov_min)/100
            cc_acov_mat(m,n) = 0;
        end
    end
end

% principal component analysis of covariance matrices (copied from the
% matlab fxn "pcacov")
[U_bb,W_bb,V_bb] = svd(bb_scov_mat); % [u,latent,pc] = svd(x);
W_bb = diag(W_bb);
[U_cc,W_cc,V_cc] = svd(cc_scov_mat); 
W_cc = diag(W_cc);

%totalvar = sum(latent);
%explained = 100*latent/totalvar;

return

%figure
%clabel(contour(cc_scov_mat,[-2000:500:2000]))
%figure
%clabel(contour(cc_acov_mat,[-2000:500:2000]))