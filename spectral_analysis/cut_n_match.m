function [spectra,baseline] = cut_n_match(spectra,aord,varargin)
% cut_n_match cuts off the baseline of both spectrum A and spectrum B and
% then resizes each spectrum so they have the same number of points

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

if ~any(strcmp(aord,{'a';'d'}))
    error('input argument aord must be either "a" or "d"');
end

% cut off baseline
baseline = cell(1,ns);
switch aord
    case 'a'
        [cut_spectra,B_bound,ind_loc] = no_baseline(spectra,aord);
        
        for s = 1:ns
            baseline{s} = spectra([1:ind_loc(s,1) ind_loc(s,2):end],s*2-1:s*2);
        end
    case 'd'
        cut_spectra = spectra(100:np-100,:);
        B_bound(1,:) = spectra(100,1:2:end);
        B_bound(2,:) = spectra(np-100,1:2:end);
        B_bound = B_bound';
        
        for s = 1:ns
            baseline{s} = spectra([1:99 np-99:end],s*2-1:s*2);
        end
end

% align cut spectra around mean magnetic field
avgB = mean(mean(B_bound,2));
B_fields = spectra(:,1:2:end);
    
for i = 1:ns
    bi = find(avgB == B_fields(:,i));
    if isempty(bi)
        k = max(find(avgB > B_fields(:,i)));
        %k = min(find(avgB < B_fields(:,i)));
        bind(i) = k; 
        B(i) = B_fields(k);
    else
        bind(i) = bi;
        B(i) = B_fields(bi);
    end
end

[np,nc] = size(cut_spectra);
dBi = bind-bind(1);
maxdi = max(dBi);
mindi = abs(min(dBi));
alined_spec = repmat(NaN,[np+maxdi+mindi nc]);
p = [1:np]';

for i = 1:ns
    bi = dBi(i);
    Ap = p-bi; 
    [c,ia,ip] = intersect(Ap,p);
    if bi < 0
        % align right
        [e,k] = setdiff(Ap,p);
        alined_spec(maxdi+ip,i*2-1:i*2)  = cut_spectra(ia,i*2-1:i*2);
        alined_spec(maxdi+ip(end)+1:maxdi+ip(end)+1+length(k)-1,i*2-1:i*2) = cut_spectra(k,i*2-1:i*2);                 
    elseif bi > 0 
        % align left
        [e,k] = setdiff(Ap,p);
        alined_spec((maxdi-length(k))+1:(maxdi-length(k))+1+length(k)-1,i*2-1:i*2) = cut_spectra(k,i*2-1:i*2);
        alined_spec(maxdi+ip,i*2-1:i*2)  = cut_spectra(ia,i*2-1:i*2);
    else
        % already aligned
        alined_spec(1+maxdi:length(ia)+maxdi,i*2-1:i*2) = cut_spectra(ia,i*2-1:i*2);
    end
end

% determine max and min magnetic fields at both high field and low field
% ends of spectra
B_bound = B_bound'; % now B_bound = 2 x ns, where top row = low field values and bottow row = high field values
maxBlow = max(B_bound(1,:));
[minBlow,sli] = min(B_bound(1,:));
[maxBhigh,shi] = max(B_bound(2,:));
minBhigh = min(B_bound(2,:));

% fill in absorbance regions that are not aligned throughout all spectra
% with zeros, this is ok since these regions are from the baseline which
% are always near zero or fluctuating around zero
B_fields = alined_spec(:,1:2:end);
bi_min_low = find(minBlow == B_fields(:,sli));
bi_max_high = find(maxBhigh == B_fields(:,shi));

for i = 1:ns
    bi_low = find(B_bound(1,i) == B_fields(:,i));
    bi_high = find(B_bound(2,i) == B_fields(:,i));
    nlow = bi_low-bi_min_low;
    nhigh = bi_max_high-bi_high;
    alined_spec(bi_min_low:bi_low-1,i*2-1:i*2) = [B_fields(bi_min_low:bi_low-1,sli) zeros(nlow,1)];
    alined_spec(bi_high+1:bi_max_high,i*2-1:i*2) = [B_fields(bi_high+1:bi_max_high,shi) zeros(nhigh,1)];
end

n = all(~isnan(alined_spec),2);
spectra = alined_spec(n,:);

return