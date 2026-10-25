function plotIntegratedFlux(x,t,phiHist,plotTitle, figure_num)

dx = x(2)-x(1);
integratedFlux = dx * sum(phiHist,1);

figure(figure_num);
semilogy(t(3:end),max(integratedFlux(3:end),realmin),'LineWidth',1.5);

xlabel('Time t');
ylabel('\int_0^L \phi(x,t) dx');
title(plotTitle);
grid on;
box on;
end