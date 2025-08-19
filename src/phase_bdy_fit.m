function output = phase_bdy_fit(S,C,V,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nrow,ncol] = size(S);

problem.objective = @fit_fxn;	
problem.x0 = [1 0.1 0.4];	 
problem.xdata = C;
ydata = S./V;
% ydata = S;
problem.ydata = ydata(:);	 
problem.lb = [0.01 C(1) C(1)];	 
problem.ub = [100 C(end) C(end)];	 
problem.solver = 'lsqcurvefit';
options = optimset('display','iter','FunValCheck','on');
problem.options = options;

% [solution,resnorm,residual,exitflag] = lsqcurvefit(problem);
[solution,fval,exitflag] = ...
    fminsearchcon(@fit_fxn,problem.x0,problem.lb,problem.ub,[],[],[],problem.options,C);

output.solution = solution;
% output.resnorm = resnorm;
% output.residual = residual;
output.fxneval = fval;
output.exitflag = exitflag;

    function F = fit_fxn(x,xdata)
        nc = length(xdata);
        B = interp1(xdata',S',[x(2) x(3)]');
        B = B';
        B = B./repmat(sum(B),nrow,1);
%         if any(isnan(B(:)))
%             keyboard;
%         end
        f = zeros(nc,2);
        for c = 1:nc
            if xdata(c) <= x(2)
                f(c,:) = [1 0];
            elseif xdata(c) > x(2) && xdata(c) < x(3)
                fb = (xdata(c)-x(2))./(x(3)-x(2));
                fa = 1-fb;
                fpb = (x(1).*fb)./(fa + x(1).*fb);
                fpa = 1-fpb;
                f(c,:) = [fpa fpb];
            elseif xdata(c) >= x(3)
                f(c,:) = [0 1];
            end
        end
%         if any(isnan(f(:)))
%             keyboard;
%         end
        Y = B*f';
%         Y = Y./V;
%         F = Y(:);
        Y = (S-Y)./V;
        F = Y(:)'*Y(:);
    end
output.B = B;
output.f = f;
end