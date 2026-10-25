function [x, t, phiHist, CHist] = runTransientSN(Nx, v, beta, lambda, sourceType, dt, g0)
%RUNTRANSIENTSN Calculate time-dependent angular flux.
%
% Nx     : Number of slabs of thickness 0.1
% v      : Unitless neutron velocity
% beta   : Delayed neutron fraction
% lambda : Delayed neutron decay constant
% sourceType    : Type of source in ['steady', 'pulse]
% dt     : Unitless time step
% g0     : Scalar incoming source strength

% Fixed material properties
Sigma_f = 0.2;
Sigma_c = 0.1;
Sigma_s = 0.1;
nu = 2;

Sigma_t = Sigma_f + Sigma_c + Sigma_s;

dx = 0.1;

% S2 quadrature
mu = [-1; 1] / sqrt(3);
w  = [1; 1];
M = length(mu);

% Spatial and temporal grids
x = ((1:Nx) - 0.5)' * dx;
Nt = round(1/dt);
t = (0:Nt) * dt;

% Initial conditions
psi = zeros(M, Nx);
phi = zeros(Nx, 1);
C = zeros(Nx, 1);

phiHist = zeros(Nx, Nt+1);
CHist = zeros(Nx, Nt+1);
phiHist(:,1) = phi;
CHist(:,1) = C;

% Iteration controls
tol = 3e-4;
maxIter = 1e5;

for n = 1:Nt

    % Boundary source for the current time step
    if strcmpi(sourceType, 'steady')
        g = g0;
    elseif strcmpi(sourceType, 'pulse')
        if n == 1
            g = g0 / dt;
        else
            g = zeros(M,1);
        end
    else
        error('sourceType must be steady or pulse.');
    end

    psiOld = psi;
    phiGuess = phi;

    for iter = 1:maxIter
        Ctrial = (C + dt * beta * nu * Sigma_f * phiGuess) / (1 + lambda * dt);
        Q = 0.5 * (Sigma_s + (1-beta)*nu*Sigma_f) * phiGuess + 0.5 * lambda * Ctrial;

        psiNew = transportSweep(psiOld, Q, mu, dx, dt, v, Sigma_t, g);
        phiNew = (w' * psiNew)';
        psiOld = psiNew;

        % Convergence check, 1e-14 to prevent divide by 0 error, finding max error in phi vector
        err = max(abs((phiNew - phiGuess) ./ max(phiNew, 1e-14)));

        if err < tol
            break
        end

        phiGuess = phiNew;
    end

    if iter == maxIter
        warning('Source iteration did not converge at t = %.6g.', t(n+1));
    end

    % Accept the converged solution
    psi = psiNew;
    phi = phiNew;
    C = Ctrial;

    phiHist(:,n+1) = phi;
    CHist(:,n+1) = C;
end
end