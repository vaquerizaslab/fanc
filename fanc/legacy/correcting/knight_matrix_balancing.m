function [x,res] = bnewt(A,tol,x0,delta,Delta,fl)
% BNEWT A balancing algorithm for symmetric matrices
%
% X = BNEWT(A) attempts to find a vector X such that
% diag(X)*A*diag(X) is close to doubly stochastic. A must
% be symmetric and nonnegative.
%
% X0: initial guess.  TOL: error tolerance.
% delta/Delta: how close/far balancing vectors can get
% to/from the edge of the positive cone.
% We use a relative measure on the size of elements.
% FL: intermediate convergence statistics on/off.
% RES: residual error, measured by norm(diag(x)*A*x - e).
% Initialise
n = size(A,1); e = ones(n,1); res=[];
if nargin < 6,  fl = 0; end
if nargin < 5,  Delta = 3; end
if nargin < 4,  delta = 0.1; end
if nargin < 3,  x0 = e; end
if nargin < 2,  tol = 1e-6; end
g=0.9; etamax = 0.1; % Parameters used in inner stopping criterion.
eta = etamax; stop_tol = tol*.5;
x = x0; rt = tol^2; v = x.*(A*x); rk = 1 - v;
rho_km1 = rk'*rk; rout = rho_km1; rold = rout;
MVP = 0;  % We'll count matrix vector products.
i = 0; % Outer iteration count.
if fl == 1, fprintf('it    in. it    res\n'), end

fprintf('starting iterations\n')

while rout > rt  % Outer iteration
    fprintf('it')
    i = i + 1; k = 0; y = e;
    innertol = max([eta^2*rout,rt]);
    while  rho_km1 > innertol %Inner iteration by CG
        k = k + 1;
        if k == 1
            Z = rk./v;  p=Z; rho_km1 = rk'*Z;
        else
            beta=rho_km1/rho_km2;
            p=Z + beta*p;
        end
        % Update search direction efficiently.
        w = x.*(A*(x.*p)) + v.*p;
        alpha = rho_km1/(p'*w);
        ap = alpha*p;

        % Test distance to boundary of cone.
        ynew = y + ap;
        if min(ynew) <= delta
            if delta == 0, break, end
            ind = find(ap < 0);
            gamma = min((delta  - y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        if max(ynew) >= Delta
            ind = find(ynew > Delta);
            gamma = min((Delta-y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        y = ynew;
        rk = rk - alpha*w;  rho_km2 = rho_km1;
        Z = rk./v;  rho_km1 = rk'*Z;
    end
    x = x.*y; v = x.*(A*x);
    rk = 1 - v; rho_km1 = rk'*rk; rout = rho_km1;
    MVP = MVP + k + 1;
    % Update inner iteration stopping criterion.
    rat = rout/rold;  rold = rout; res_norm = sqrt(rout);
    eta_o = eta;  eta = g*rat;
    if g*eta_o^2 > 0.1
        eta = max([eta,g*eta_o^2]);
    end
    eta = max([min([eta,etamax]),stop_tol/res_norm]);
    if fl == 1
        fprintf('%3d %6d   %.3e %.3e %.3e \n', i,k, r_norm,min(y),min(x));
        res=[res; r_norm];
    end
end
fprintf('Matrix-vector products = %6d\n', MVP)
end