function [B,I,indices] = peak(spectra,B0)
% find magnetic field value of nearest peak to B0 field value

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

for i =1:ns
    B_field = spectra(:,i*2-1);
    spectrum = spectra(:,i*2);
    
    % take an initial step in the direction of the gradient
    ind0 = find(B0 == B_field);
    if isempty(ind0)
        j = max(find(B0 > B_field));
        k = min(find(B0 < B_field));
        ind0 = j; % or k i guess
    end
    
    G = gradient(spectrum);
    GatB0 = G(ind0);
    top = false;
    switch sign(GatB0)
        case 1 % GatB0 > 0
            % intensity increases with increasing field
            indb = ind0;
            Gb = G(indb);
            inda = indb+1;
            Ga = G(inda);
        case 0 % GatB0 == 0
            % already at nearest peak
            Gb = 0;
            Ga = 0;
            top = true;
        case -1 % GatB0 < 0
            % intensity decreases with increasing field
            indb = ind0;
            Gb = G(indb);
            inda = indb-1;
            Ga = G(inda);
        otherwise
            error(sprintf('The %d-th spectrum has an undefined gradient at B-field index %d',i,ind0));
    end
    
    if sign(Gb) ~= sign(Ga)
        % find index that has a higher value 
        if spectrum(inda) > spectrum(indb)
            indf = inda;
        else
            indf = indb;
        end
        
        top = true;
    end
    
    % find top of nearest peak following the gradient
    while ~top % (sign(Gb) == sign(Ga)) == (Gb < 0 & Ga > 0) | (Gb > 0 & Ga < 0)
        if all(sign([Gb Ga]) == 1)
            indb = inda;
            Gb = G(indb);
            inda = indb+1;
            Ga = G(inda);
        elseif all(sign([Gb Ga]) == -1)
            indb = inda;
            Gb = G(indb);
            inda = indb-1;
            Ga = G(inda);
        elseif all(sign([Gb Ga]) == 0)
            indf = inda;
            top = true;
        elseif any([Gb Ga] == 0)
            if Gb == 0
                indf = indb;
            else
                indf = inda;
            end
            
            top = true;
        elseif sign(Gb) ~= sign(Ga)
            % find index that has a higher value 
            if spectrum(inda) > spectrum(indb)
                indf = inda;
            else
                indf = indb;
            end
            
            top = true;
        end
    end
    
    B(i) = B_field(indf);
    I(i) = spectrum(indf);
    indices(i) = indf;
end

return