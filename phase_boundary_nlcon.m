function [nlin,nleq,Dnlin,Dnleq] = phase_boundary_nlcon(x,S,C,e_spec)
% x = [Kp C_alpha C_beta]
% S = spectra matrix
% C = spectal composition vector
% e_spec = experimental or data spectra vector (each spectrum stacked on
% top of one another)

[nb,ns] = size(S);
nc = ns;

nlin = [];
%nlin = [-((x(1).*(x(3)-C'))./(C'-x(2)+(x(1).*(x(3)-C'))));...
%        ((x(1).*(x(3)-C'))./(C'-x(2)+(x(1).*(x(3)-C'))))-1;...
%        -((C'-x(2))./(C'-x(2)+(x(1).*(x(3)-C'))));...
%        ((C'-x(2))./(C'-x(2)+(x(1).*(x(3)-C'))))-1];          % Nonlinear inequalities at x, nlin <= 0
nleq = ((x(1).*(x(3)-C'))./(C'-x(2)+(x(1).*(x(3)-C'))))+((C'-x(2))./(C'-x(2)+(x(1).*(x(3)-C'))))-1;  % Nonlinear equalities at x, nleq = 0

if nargout > 2   % nonlcon called with 4 outputs
   %Dnlin = zeros(3,4*nc);      % Gradients of the inequalities
   Dnleq = zeros(3,nc);    % Gradients of the equalities
   % calculate gradient D
   kp = x(1);
   a = x(2);
   b = x(3);
   c = C;
   d1kp = (((b-c)./(c-a+kp.*(b-c)))-(kp.*(b-c).^2./(c-a+kp.*(b-c)).^2));
   d1a = (kp.*(b-c)./(c-a+kp.*(b-c)).^2);
   d1b = ((kp./(c-a+kp.*(b-c)))-(kp.^2.*(b-c)./(c-a+kp.*(b-c)).^2));
   d2kp = -((c-a).*(b-c)./(c-a+kp.*(b-c)).^2);
   d2a = (-(1./(c-a+kp.*(b-c)))+((c-a)./(c-a+kp.*(b-c)).^2));
   d2b = -((c-a).*kp./(c-a+kp.*(b-c)).^2);
   %Dnlin = [-d1kp,d1kp,-d2kp,d2kp;...             
   %      -d1a,d1a,-d2a,d2a;...      
   %      -d1b,d1b,-d2b,d2b];
   Dnlin = [];
   Dnleq = [(d1kp+d2kp);(d1a+d2a);(d1b+d2b)];
end