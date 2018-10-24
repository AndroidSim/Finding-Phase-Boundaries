function spectra = align_spectra(spectra)
% align_spectra: aligns spectra 
%   spectra = (number magnetic fields x 2*number of spectra)
% syntac: alined_spec = align_spectra(spectra)

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
B_fields = spectra(:,1:2:end);
rB_fields = round(B_fields*100)/100;

for s = 1:ns-1
    if isequal(rB_fields(:,s),rB_fields(:,s+1))
        alined(s) = true;
    else
        alined(s) = false;
    end
end

if all(alined)
    % all spectra are already aligned
    return;
end
    
lowB = B_fields(1,:);
highB = B_fields(end,:);
ranges = highB-lowB;
dranges = max(ranges)-ranges;
Binc = ranges/np;

% 2 ways to align spectra:
%   1)  if the B field values and increments not the same for all spectra then must
%   align them by interpolation of a constant range and increment.
%   2)  if the B field values and increments the same then they can be
%   aligned by shifting the spectra (faster than interpolation).

if ~all(diff(Binc) == 0)
    minlowB = min(lowB);
    maxhighB = max(highB);
    avg_inc = mean(Binc);
    Bi = [minlowB:avg_inc:maxhighB]';
    S = spectra(:,2:2:end);
    
    for s = 1:ns
        ispectra(:,s*2-1) = Bi;
        ispectra(:,s*2) = interp1(B_fields(:,s),S(:,s),Bi,'linear',0);
    end

    spectra = ispectra; 
else
    [np,nc] = size(spectra);
    medB = median(median(B_fields));
    medB = round(medB);

    for s = 1:ns
        Bi = find(medB == B_fields(:,s));
        if isempty(Bi)
            k = max(find(medB > B_fields(:,s)));
            %k = min(find(avgB < B_fields(:,i)));
            Bind(s) = k; 
            B(s) = B_fields(k,s);
        else
            Bind(s) = Bi;
            B(s) = B_fields(Bi,s);
        end
    end

    dBi = Bind-Bind(1);
    maxdi = max(dBi);
    mindi = abs(min(dBi));
    alined_spec = repmat(NaN,[maxdi+np+mindi nc]);

    for s = 1:ns
        bi = dBi(s);
        alined_spec(maxdi-bi+1:end-mindi-bi,s*2-1:s*2) = spectra(:,s*2-1:s*2);
    end

    B_fields = alined_spec(:,1:2:end);
    B_bound(1,:) = lowB;
    B_bound(2,:) = highB;

    % determine max and min magnetic fields at both high field and low field
    % ends of spectra
    maxBlow = max(B_bound(1,:));
    [minBlow,sli] = min(B_bound(1,:));
    [maxBhigh,shi] = max(B_bound(2,:));
    minBhigh = min(B_bound(2,:));

    % fill in absorbance regions that are not aligned throughout all spectra
    % with zeros, this is ok since these regions are from the baseline which
    % are always near zero or fluctuating around zero
    bi_min_low = find(minBlow == B_fields(:,sli));
    bi_max_high = find(maxBhigh == B_fields(:,shi));

    for s = 1:ns
        bi_low = find(B_bound(1,s) == B_fields(:,s));
        bi_high = find(B_bound(2,s) == B_fields(:,s));
        nlow = bi_low-bi_min_low;
        nhigh = bi_max_high-bi_high;
        alined_spec(bi_min_low:bi_low-1,s*2-1:s*2) = [B_fields(bi_min_low:bi_low-1,sli) zeros(nlow,1)];
        alined_spec(bi_high+1:bi_max_high,s*2-1:s*2) = [B_fields(bi_high+1:bi_max_high,shi) zeros(nhigh,1)];
    end

    spectra = alined_spec;
end

return