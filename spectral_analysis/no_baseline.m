function [cut_specs,B_bound,ind_loc] = no_baseline(spectra,aord)
% no_baseline cuts off the baseline and leaves just the peaks and
% valleys of the spectrum
%
% aord = 'a' if absorbance spectra, or 'd' for derivative

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field intensity_values]');
end

ns = nc/2;
% Bc = mean(mean(spectra(:,1:2:end)));
Bc = 3326;

if ~any(strcmp(aord,{'a';'d'}))
    error('input argument aord must be either "a" or "d"');
end

cut_specs = repmat(NaN,size(spectra));

if aord == 'a'
    for i = 1:ns
        B_field = spectra(:,i*2-1);
        spectrum = spectra(:,i*2);
        
        % find central peak
        [B,index] = peak([B_field spectrum],Bc);
        
        % find neighboring valleys
        stepleft = index-10;
        stepright = index+10;
        Bleft = B_field(stepleft);
        Bright = B_field(stepright);
        [Bleft,ileft] = valley([B_field spectrum],Bleft);
        [Bright,iright] = valley([B_field spectrum],Bright);
        
        % find neighboring peaks to above valleys
        stepleft = ileft-10;
        stepright = iright+10;
        Bleft = B_field(stepleft);
        Bright = B_field(stepright);
        [Bleft,ileft] = peak([B_field spectrum],Bleft);
        [Bright,iright] = peak([B_field spectrum],Bright);
        
        % finding the valley from above peaks gives the boundary B_field
        % values, stored in B_bound
        stepleft = ileft-10;
        stepright = iright+10;
        Bleft = B_field(stepleft);
        Bright = B_field(stepright);
        
        ileftb = stepleft;
        ilefta = stepleft-1;
        vb = spectrum(ileftb);
        va = spectrum(ilefta);
        
        while sign(va) == sign(vb)
            ileftb = ileftb-1;
            ilefta = ilefta-1;
            
            if ilefta == 1
                break
            else
                vb = spectrum(ileftb);
                va = spectrum(ilefta);
            end
        end
        
        Bleft = B_field(ileftb);
        
        irightb = stepright;
        irighta = stepright+1;
        vb = spectrum(irightb);
        va = spectrum(irighta);
        
        while sign(va) == sign(vb)
            irightb = irightb+1;
            irighta = irighta+1;
            
            if irighta == np
                break
            else
                vb = spectrum(irightb);
                va = spectrum(irighta);
            end
        end
        
        Bright = B_field(irightb);
        
        %[Bleft,ileft] = valley([B_field spectrum],Bleft);
        %[Bright,iright] = valley([B_field spectrum],Bright);
        
        indices = find(Bleft <= B_field & B_field <= Bright);
        cut_specs(indices,i*2-1:i*2) = [B_field(indices) spectrum(indices)];
        B_bound(i,:) = [Bleft Bright]; 
        ind_loc(i,:) = [indices(1) indices(end)];
    end
end

if aord == 'd'
    for i = 1:ns
        B_field = spectra(:,i*2-1);
        spectrum = spectra(:,i*2);
        
        % find central peak
        [B,index] = peak([B_field spectrum],Bc);
        
        % find neighboring valleys
        stepleft = index-10;
        stepright = index+10;
        Bleft = B_field(stepleft);
        Bright = B_field(stepright);
        [Bleft,ileft] = valley([B_field spectrum],Bleft);
        [Bright,iright] = valley([B_field spectrum],Bright);
        
        % find neighboring peaks to above valleys
        stepleft = ileft-10;
        stepright = iright+10;
        Bleft = B_field(stepleft);
        Bright = B_field(stepright);
        [Bleft,ileft] = peak([B_field spectrum],Bleft);
        [Bright,iright] = peak([B_field spectrum],Bright);
        
        % finding the valley on the left gives the left B_field boundary
        stepleft = ileft-10;
        Bleft = B_field(stepleft);
        %[Bleft,ileft] = valley([B_field spectrum],Bleft);
        
        ileftb = stepleft;
        ilefta = stepleft-1;
        vb = spectrum(ileftb);
        va = spectrum(ilefta);
        
        while sign(va) == sign(vb)
            ileftb = ileftb-1;
            ilefta = ilefta-1;
            
            if ilefta == 1
                break
            else
                vb = spectrum(ileftb);
                va = spectrum(ilefta);
            end
        end
        
        Bleft = B_field(ileftb);
        
        stepright = iright+10;
        Bright = B_field(stepright);
        [Bright,iright] = valley([B_field spectrum],Bright);
        
        % and finding the peak on the right gives the right B_field
        % boundary
        stepright = iright+10;
        Bright = B_field(stepright);
        %[Bright,iright] = peak([B_field spectrum],Bright);
        
        irightb = stepright;
        irighta = stepright+1;
        vb = spectrum(irightb);
        va = spectrum(irighta);
        
        while sign(va) == sign(vb)
            irightb = irightb+1;
            irighta = irighta+1;
            
            if irighta == np
                break
            else
                vb = spectrum(irightb);
                va = spectrum(irighta);
            end
        end
        
        Bright = B_field(irightb);
        
        indices = find(Bleft <= B_field & B_field <= Bright);
        cut_specs(indices,i*2-1:i*2) = [B_field(indices) spectrum(indices)];
        B_bound(i,:) = [Bleft Bright];
        ind_loc(i,:) = [indices(1) indices(end)];
    end
end

% find(any(~isnan(ans'),2))
return