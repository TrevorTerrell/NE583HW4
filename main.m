%% Prep
clear; clc; close all;

dx = 0.1;
dt = 1/6000;

% Use the number of slabs obtained in Assignment 2.
Nx_A2 = 30;

% S2 incoming source convention from Assignment 2
mu = [-1; 1] / sqrt(3);
g0 = zeros(2,1);
g0(mu > 0) = 1 ./ mu(mu > 0);

%% Question 1: Steady source, prompt neutrons only
[x1,t1,phi1,~] = runTransientSN(Nx_A2, 20, 0, 0, 'steady', dt, g0);

plotFluxSnapshots(x1,t1,phi1, 'Q1: Steady source, no delayed neutrons', 1);

%% Question 2: Pulsed source, prompt neutrons only
[x2,t2,phi2,~] = runTransientSN( Nx_A2, 20, 0, 0, 'pulse', dt, g0);

plotFluxSnapshots(x2,t2,phi2, 'Q2: Pulsed source, no delayed neutrons', 2);

%% Question 3: Pulsed source with one delayed group
beta = 0.0065;
lambda = 1;

[x3,t3,phi3,C3] = runTransientSN(64, 1000, beta, lambda, 'pulse', dt, g0);

plotFluxSnapshots(x3,t3,phi3, 'Q3: Pulsed source with delayed neutrons', 3);

plotIntegratedFlux(x3,t3,phi3, 'Q3: Integrated scalar flux', 4);

%% Question 4: 65 slabs, sweep over beta
Nx4 = 65;
betaValues = linspace(0.001,0.01,7);

figure(5);

for k = 1:length(betaValues)

    beta = betaValues(k);

    [x4,t4,phi4,~] = runTransientSN(Nx4, 1000, beta, 1, 'pulse', dt, g0);

    integratedFlux = dx * sum(phi4,1);

    semilogy(t4(3:end), max(integratedFlux(3:end), realmin), 'LineWidth', 1.4, 'DisplayName', sprintf('\\beta = %.4f',beta));
    hold on;

    % % Show the spatial-time evolution for one representative beta.
    % if abs(beta - 0.0065) < 1e-10
    %     figure(5);
    %     surf(x4,t4,phi4');
    %     shading interp;
    %     xlabel('Position x');
    %     ylabel('Time t');
    %     zlabel('Scalar flux \phi(x,t)');
    %     title('Q4: Scalar flux evolution, \beta = 0.0065');
    %     colorbar;
    %     view(45,30);
    % end
end

figure(5)%findobj('Type','figure','-and','Number',1));
% If figure selection is inconvenient, create the integrated-flux
% figure separately and run the beta loop there.

xlabel('Time t');
ylabel('Spatially integrated scalar flux');
title('Q4: Integrated flux for different delayed fractions');
legend('Location','best');
grid on;
box on;