function [best_p,best_chisq] = boundary_fit(xdata,ydata,bdy,varargin)

% lb = [0;0;0;0;.5];
% ub = [1;1;1;1;2];
% options = optimset('display','iter','maxfunevals',2000);
% [p,chisq] = lsqcurvefit(@boundary_fit_fxn,[0.1 0.5 0.4 0.7 1]',s_xdata,s_ydata,lb,ub,options,bdy);

best_p = zeros(5,5);
best_chisq = zeros(5,1);
nfilled = 1;

if all(size(bdy) > 1) & ndims(bdy) == 2 % ie if matrix
    [nbp,nbc] = size(bdy); % nbp = number of boundary pts, nc = number of components
    ncv = nbc-1; % ncv = number of compositional variables
    % make boundary convex
    bdy = make_convex_bdy(bdy);
    % declare compositional basis vectors
    % for ternary phase diagram using Xsm and Xchol
    es = [1 0]; % [x y]
    ec = [cos(pi/3) sin(pi/3)]; % [x y]
    bdy = bdy(:,1)*es+bdy(:,3)*ec;
    b = linspace(0,1,2*nbp)'; % b = [0:0.02:1]';
    [bdy_length,bdy_intervals] = bdy_fxn(bdy,'linear');
    bdy = interp1(bdy_intervals(:,1),bdy,b,'pchip');
else
    error('1st argument must be an nx3 matrix (Xsm Xdopc Xchol) specifying boundary points');
end

traj = varargin{1};
bt = bdypt2b(linspace(0,1,length(bdy))',bdy,tern2cart(traj,1))
%traja = [0:0.05:1];
%trajb = [0:0.05:1];
cptc = [0:0.05:1];
cptd = [0:0.05:1];
kp = [0.5:0.1:0.9 1:0.5:3];

bdy = cart2tern(bdy,1);

%for a = 1:length(traja)
    %for b = 1:length(trajb)
        for c = 1:length(cptc)
            for d = 1:length(cptd)
                for p = 1:length(kp)
                    %if traja(a) == trajb(b) | (any(traja(a) == [0 1]) & any(trajb(b) == [0 1]))
                    %    continue;
                    %end
                    
                    if cptc(c) == cptd(d) | (any(cptc(c) == [0 1]) & any(cptd(d) == [0 1]))
                        continue;
                    end
                    
                    %ytheor = feval(@boundary_fit_fxn,[traja(a);trajb(b);cptc(c);cptd(d);kp(p)],xdata,bdy);
                    disp([cptc(c) cptd(d) kp(p)]);
                    ytheor = feval(@boundary_fit_fxn,[bt;cptc(c);cptd(d);kp(p)],xdata,bdy);
                    chisq = norm(ytheor-ydata);
                    params = [bt;cptc(c);cptd(d);kp(p)];
                    
                    % update stack of best fit parameters
                    if (isempty(find(best_chisq == 0)))
                        [best_chisq,index] = sort(best_chisq);
                        best_p = best_p(index,:);
                        i = find(chisq <= best_chisq);
                        if (~isempty(i))
                            i = max(i);
                            best_chisq(i) = chisq;
                            %best_p(i,:) = [traja(a) trajb(b) cptc(c) cptd(d) kp(p)];
                            best_p(i,:) = params';
                        end
                    else
                        best_chisq(nfilled) = chisq;
                        %best_p(nfilled,:) = [traja(a) trajb(b) cptc(c) cptd(d) kp(p)];
                        best_p(nfilled,:) = params';
                        nfilled = nfilled+1;
                    end
                end
            end
        end
        %end
    %end

return