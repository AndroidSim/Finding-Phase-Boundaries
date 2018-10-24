function met_mat = metric_matrix(spectra,aord,met_type)
% metric_matric calculates a matrix where every element on and below the
% main diagonal contains the metric between two spectra
%
% spectra = matrix of the spectra to compare
% aord = 'a' for absorbance spectra or 'd' for derivative spectra
% met_type = metric type, see metric.m or pdist.m

[np,nc] = size(spectra);

%if nc < 2 || rem(nc,2) ~= 0
%    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
%end

ns = np; %nc/2;
met_mat = zeros([ns ns]);
%spectra = cut_n_match(spectra,aord);
%spectra = spectra(:,2:2:end)';

for i = 1:ns
    if isequal(i,1)
        spectra(i,:) = mean(spectra(i:i+2,:));
    elseif isequal(i,ns)
        spectra(i,:) = mean(spectra(i-2:i,:));
    else
        spectra(i,:) = mean(spectra(i-1:i+1,:));
    end
end

met_mat = squareform(pdist(spectra,'correlation'));

%for i = 1:ns-1
%    for j = 1:i
%        if isequal(i,1) & ~isequal(j,1)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i:i+2,:),spectra(j-1:j+1,:));
%        elseif ~isequal(i,1) & isequal(j,1)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i-1:i+1,:),spectra(j:j+2,:));
%        elseif isequal(i,1) & isequal(j,1)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i:i+2,:),spectra(j:j+2,:));
%        elseif isequal(i,ns) & ~isequal(j,ns)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i-2:i,:),spectra(j-1:j+1,:));
%        elseif ~isequal(i,ns) & isequal(j,ns)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i-1:i+1,:),spectra(j-2:j,:));
%        elseif isequal(i,ns) & isequal(j,ns)
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i-2:i,:),spectra(j-2:j,:));
%        else
%            met_mat(i,j) = feval(@my_metric_fun,spectra(i-1:i+1,:),spectra(j-1:j+1,:));
%        end
%    end
%end

%for i = 1:ns-1
%    for j = 1:i
%        met_mat(i,j) = feval(@my_metric_fun,spectra(i,:),spectra(j,:));
%    end
%end

%met_mat = met_mat+met_mat';
% use diag() here?
%for i = 1:ns
%    met_mat(i,i) = met_mat(i,i)/2;
%end

return