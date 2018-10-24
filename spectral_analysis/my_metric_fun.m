function d = my_metric_fun(a,b)
% this test metric is a variant of the correlation metric of the 'pdist'
% matlab function

[ra,ca] = size(a);
[rb,cb] = size(b);

if ~isequal(ca,cb)
    error('the number of columns of the data matrices must be equal');
end

xrbar = mean(a);
xsbar = mean(b);
xr = a(2,:);
xs = b(2,:);
xrcorr = xr - xrbar;
xscorr = xs - xsbar;

d = 1 - ((xrcorr*xscorr')/(sqrt(xrcorr*xrcorr')*sqrt(xscorr*xscorr')));