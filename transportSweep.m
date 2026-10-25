function psiNew = transportSweep(psiOld, Q, mu, dx, dt, v, Sigma_t, g)
%TRANSPORTSWEEP Perform one implicit S2 transport sweep.
%
% psiOld : M-by-Nx angular flux at previous time
% Q      : Nx-by-1 isotropic internal source
% mu     : M-by-1 angular ordinates
% g      : M-by-1 incoming boundary angular flux
%
% Vacuum boundary conditions apply where g = 0.

[M, Nx] = size(psiOld);
psiNew = zeros(M, Nx);

timeCoeff = dx / (v * dt);

for m = 1:M

    am = abs(mu(m));
    denom = timeCoeff + am + Sigma_t * dx;

    if mu(m) > 0
        % Sweep from left to right
        psiIn = g(m);

        for i = 1:Nx
            psiNew(m,i) = (timeCoeff * psiOld(m,i) + am * psiIn + dx * Q(i)) / denom;

            psiIn = psiNew(m,i);
        end

    else
        % Sweep from right to left
        psiIn = g(m);

        for i = Nx:-1:1
            psiNew(m,i) = (timeCoeff * psiOld(m,i) + am * psiIn + dx * Q(i)) / denom;

            psiIn = psiNew(m,i);
        end
    end
end
end