function [Bpeaks,Bvalleys,Ipeaks,Ivalleys,Vpeaks,Vvalleys] = get_peaks_valleys(spectra,aord)
% no_baseline cuts off the baseline and leaves just the peaks and
% valleys of the spectrum
%
% aord = 'a' if absorbance spectra, or 'd' for derivative

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field intensity_values]');
end

ns = nc/2;
Bc = mean(mean(spectra(:,1:2:end)));

switch aord
    case 'a'
        Bpeaks = repmat(NaN,[3 ns]);
        Bvalleys = repmat(NaN,[2 ns]);
        Ipeaks = repmat(NaN,[3 ns]);
        Ivalleys = repmat(NaN,[2 ns]);
        Vpeaks = repmat(NaN,[3 ns]);
        Vvalleys = repmat(NaN,[2 ns]);
        
        for i = 1:ns
            B_field = spectra(:,i*2-1);
            spectrum = spectra(:,i*2);
        
            % find central peak
            [B,index] = peak([B_field spectrum],Bc);
            Bpeaks(2,i) = B;
            Ipeaks(2,i) = index;
            Vpeaks(2,i) = spectrum(index);
        
            % find neighboring valleys
            stepleft = index-10;
            stepright = index+10;
            Bleft = B_field(stepleft);
            Bright = B_field(stepright);
            [Bleft,ileft] = valley([B_field spectrum],Bleft);
            [Bright,iright] = valley([B_field spectrum],Bright);
            Bvalleys(:,i) = [Bleft; Bright];
            Ivalleys(:,i) = [ileft; iright];
            Vvalleys(:,i) = spectrum(Ivalleys(:,i));
        
            % find neighboring peaks to above valleys
            stepleft = ileft-10;
            stepright = iright+10;
            Bleft = B_field(stepleft);
            Bright = B_field(stepright);
            [Bleft,ileft] = peak([B_field spectrum],Bleft);
            [Bright,iright] = peak([B_field spectrum],Bright);
            Bpeaks(1,i) = Bleft;
            Bpeaks(3,i) = Bright;
            Ipeaks(1,i) = ileft;
            Ipeaks(3,i) = iright;
            Vpeaks(1,i) = spectrum(ileft);
            Vpeaks(3,i) = spectrum(iright);
        end
    case 'd'
        Bpeaks = repmat(NaN,[3 ns]);
        Bvalleys = repmat(NaN,[3 ns]);
        Ipeaks = repmat(NaN,[3 ns]);
        Ivalleys = repmat(NaN,[3 ns]);
        
        for i = 1:ns
            B_field = spectra(:,i*2-1);
            spectrum = spectra(:,i*2);
        
            % find central peak
            [B,index] = peak([B_field spectrum],Bc);
            Bpeaks(2,i) = B;
            Ipeaks(2,i) = B;
        
            % find neighboring valleys
            stepleft = index-5;
            stepright = index+5;
            Bleft = B_field(stepleft);
            Bright = B_field(stepright);
            [Bleft,ileft] = valley([B_field spectrum],Bleft);
            [Bright,iright] = valley([B_field spectrum],Bright);
            Bvalleys(1:2,i) = [Bleft; Bright];
            Ivalleys(1:2,i) = [ileft; iright];
        
            % find neighboring peaks to above valleys
            stepleft = ileft-5;
            stepright = iright+5;
            Bleft = B_field(stepleft);
            Bright = B_field(stepright);
            [Bleft,ileft] = peak([B_field spectrum],Bleft);
            [Bright,iright] = peak([B_field spectrum],Bright);
            Bpeaks(1,i) = Bleft;
            Bpeaks(3,i) = Bright;
            Ipeaks(1,i) = ileft;
            Ipeaks(3,i) = iright;
        
            stepright = iright+5;
            Bright = B_field(stepright);
            [Bright,iright] = valley([B_field spectrum],Bright);
            Bvalleys(3,i) = Bright;
            Ivalleys(3,i) = iright;
        end
    otherwise
        error('input argument aord must be either "a" or "d"');
end

% find(any(~isnan(ans'),2))
return