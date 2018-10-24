function ohs = outer_hyperfine_splittings(spectra)
% call syntax: ohs = outer_hyperfine_splittings(spectra)
% calculates the outer hyperfine splitting for the spectra after asking for
% the low field and high field ranges

[np,nc] = size(spectra);
if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field signal_values]');
end
ns = nc/2;

B = spectra(:,1:2:end);
I = spectra(:,2:2:end);

for s = 1:ns
    disp(sprintf('click 2 points on the low field side to indicate the possible range:\n'));
    plot_spectra([B(:,s) I(:,s)]);
    [x,y] = ginput(2); % expect x = magnetic field, y = intensity
    pause(2);
    [c,ix1,iB1] = intersect(round(x(1)*10)/10,round(B(:,1).*10)./10); % intersect(x,B(:,s));
    [c,ix2,iB2] = intersect(round(x(2)*10)/10,round(B(:,1).*10)./10);
    [maxIlow,imaxI] = max(I(iB1:iB2,s)); 
    BmaxIlow = B(iB1+imaxI-1,s);
    disp(sprintf('click 2 points on the high field side to indicate the possible range:\n'));
    [x,y] = ginput(2); % expect x = magnetic field, y = intensity
    pause(2);
    close;
    [c,ix1,iB1] = intersect(round(x(1)*10)/10,round(B(:,1).*10)./10); % intersect(x,B(:,s));
    [c,ix2,iB2] = intersect(round(x(2)*10)/10,round(B(:,1).*10)./10);
    [minIhigh,iminI] = min(I(iB1:iB2,s)); 
    BminIhigh = B(iB1+iminI-1,s);
    ohs(s,1) = BminIhigh-BmaxIlow;
end

return