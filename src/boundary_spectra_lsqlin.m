function [D,D_con,Sbar,x,A] = boundary_spectra_lsqlin(spectra,C)
% spectra should be derivative spectra from .dat files converted from ascii
% (.asc) files

[nb,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field intensity_values]');
end

ns = nc/2;

% convert derivative spectra to absorbance spectra
spectra = deriv2abs(spectra);

% normalize absorbance spectra
spectra = normalize_spectra(spectra,'a');

% cut off baseline and match length of cut spectra
[spectra,baseline] = cut_n_match(spectra,'a');

% get baseline statistics
avg_bline = repmat(NaN,ns,1);
std_bline = repmat(NaN,ns,1);
highest_pos = repmat(NaN,ns,1);
lowest_neg = repmat(NaN,ns,1);

for s = 1:ns
    if isempty(baseline{s}), continue, end
    avg_bline(s) = mean(baseline{s}(:,2));
    std_bline(s) = std(baseline{s}(:,2),1);
    
    temp = max(baseline{s}(baseline{s}(:,2) > 0,2));
    if ~isempty(temp)
        highest_pos(s) = temp;
    end
    
    temp = min(baseline{s}(baseline{s}(:,2) < 0,2));
    if ~isempty(temp)
        lowest_neg(s) = temp;
    end
end

% interpolate spectra for constant magnetic field increments
interval = (spectra(end,1) - spectra(1,1))/length(spectra(:,1));
interval = round(interval*1000)/1000;
Bi = [spectra(1,1):interval:spectra(end,1)]';
ispectra = zeros([length(Bi) nc]);

for s = 1:ns
    B = spectra(:,s*2-1);
    
    if any(diff(B) == 0)
        i = find(diff(B) == 0);
        
        for k = 1:length(i)
            inc = (B(i(k)+2)-B(i(k)-1))/3;
            B(i(k)) = B(i(k)-1)+inc;
            B(i(k)+1) = B(i(k))+inc;
        end
    end
    
    if any(diff(B) < 0)
        i = find(diff(B) < 0);
        
        for k = 1:length(i)
            low = 1;
            high = 2;
            m = (B(i(k)+high)-B(i(k)-low))/(high+low);
            
            while m <= 0
                low = low-1;
                high = high+1;
                m = (B(i(k)+high)-B(i(k)-low))/(high+low);
            end
            
            B(i(k)) = B(i(k)-1)+m;
            B(i(k)+1) = B(i(k))+m;
        end
    end
    
    A = spectra(:,s*2);
    Ai = interp1(B,A,Bi); % default is linear
    ispectra(:,s*2-1:s*2) = [Bi Ai];
end

spectra = ispectra;

% normalize spectra so that the sum of absorbances equals one
for s = 1:ns
    spectra(:,s*2) = spectra(:,s*2)./sum(spectra(:,s*2));
end

% reduce spectra to just absorbance values and subtract mean spectrum
S = spectra(:,2:2:end);
% S = S - repmat(mean(S,2),1,ns);

% perform singular value decomposition to obtain eigenspectra basis
[U,W,V] = svd(S,0);
B = U(:,1:2); % B = 2 component eigenbasis

% solve for eigenspectra coefficient matrix D by linear least squares
options = optimset('display','off','largescale','off');
warning off;

% constraints
Aeq = [sum(U(:,1)) sum(U(:,2))];
beq = 1;
A = -B;
b = zeros(size(U(:,1)));
D = W(1:2,1:2)*V(:,1:2)';

for s = 1:ns
    % [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq,[],[],x0(:,s)); % ...[],[],x0,options
    [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq);
    % x = lsqlin(B,S(:,s),[],[],Aeq,beq);
    % x = lsqlin(B,S(:,s),A,b);
    % x = B\S(:,s);
    D_con(:,s) = x;
end
 
Sbar = B*D_con;
D_con = D_con';
D = D';

% solve for boundary coordinates a_alpha,b_alpha,a_beta,and b_beta using D
% D = [a b]; x = [kp a_alpha a_beta c_alpha c_beta];

options = optimset('display','iter','largescale','on');
maxa = max(D_con(:,1));
mina = min(D_con(:,1));
[x,resnorm] = lsqcurvefit(@boundary_spectra_fun,[1 D_con(1,1) D_con(end,1) C(1) C(end)],C,D_con(:,1),[0.1 mina mina 0 0],[10 maxa maxa 1 1],options);
A = feval(@boundary_spectra_fun,x,C);

figure;
plot([D_con(:,1) A]);

% constraints
%Aeq = [sum(U(:,1)) sum(U(:,2)) 0 0;0 0 sum(U(:,1)) sum(U(:,2))];
%beq = [1;1];
%A = [1 0];
%b = zeros(size(U(:,1)));

%options = optimset('display','iter','largescale','on');
%x0 = [D(9,1) D(9,2) D(10,1) D(10,2)]';
%[x,fval] = fsolve(@boundary_spectra_fun,x0,options,D(:,1),D(:,2));
%[x,fval] = fmincon(@boundary_spectra_fun,x0,[],[],Aeq,beq,[],[],[],options,D(:,1),D(:,2));

warning on;

return