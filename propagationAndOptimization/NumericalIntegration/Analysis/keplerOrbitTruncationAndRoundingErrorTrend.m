set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultTextInterpreter','latex'); 

close all
clear all
clc

folder = '../../SimulationOutput/IntergratorAndPropagatorInfluence/';
saveResults = true;

load(strcat(folder,'integrationErrorBehaviour.dat'));

subplot(1,2,1)
loglog(10.^integrationErrorBehaviour(1:45,1),sqrt(sum(integrationErrorBehaviour(1:45,2:4).^2')),'LineWidth',2)
grid on
xlim([1E-2,1E2])
xlabel('Time step [s]')
ylabel('Position error norm after 2.5 hours[m]')

subplot(1,2,2)
loglog(10.^integrationErrorBehaviour(1:45,1),sqrt(sum(integrationErrorBehaviour(1:45,5:7).^2')),'LineWidth',2)
grid on
xlim([1E-2,1E2])
xlabel('Time step [s]')
ylabel('Velocity error norm after hours[m]')

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

if( saveResults )
    saveas(gcf,'truncationRoundingErrorRk4','png');
end

