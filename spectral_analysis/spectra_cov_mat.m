function cov_mat = spectra_cov_mat(spectra)
% spectra must be aligned and only contain the intensity values

spectra = spectra';
[ns,np] = size(spectra);
spectra = diff(spectra);
 
% calculate covariance matrix
cov_mat = spectra' * spectra / (ns-1);
return