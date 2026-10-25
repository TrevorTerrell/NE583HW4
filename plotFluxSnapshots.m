function plotFluxSnapshots(x,t,phiHist,plotTitle, figure_num)

timeTargets = linspace(0,1,7);
timeIndices = zeros(size(timeTargets));

for k = 1:length(timeTargets)
    [~,timeIndices(k)] = min(abs(t-timeTargets(k)));
end

figure(figure_num);
hold on;

for k = 1:length(timeIndices)
    j = timeIndices(k);

    plot(x,phiHist(:,j),'LineWidth',1.5, ...
        'DisplayName',sprintf('t = %.4f',t(j)));
end

xlabel('Position x');
ylabel('Scalar flux \phi(x,t)');
title(plotTitle);
legend('Location','best');
grid on;
box on;
hold off;
end