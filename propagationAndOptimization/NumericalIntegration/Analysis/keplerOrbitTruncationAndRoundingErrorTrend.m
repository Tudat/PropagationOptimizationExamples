set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultTextInterpreter','latex'); 

set(0,'defaultAxesFontSize',14)
set(0,'defaultTextFontSize',14)

close all
clear all
clc

folder = '../../SimulationOutput/NumericalIntegration/';
saveResults = true;
saveFolder = '/home/dominic/Documents/Courses/Propagation and Optimization/IntegrationFigures/';

load(strcat(folder,'integrationErrorBehaviour.dat'));
load(strcat(folder,'integrationErrorBehaviour_long.dat'));

loglog(10.^integrationErrorBehaviour(:,1),sqrt(sum(integrationErrorBehaviour(:,2:4).^2')),'-o','LineWidth',2)
grid on
xlim([1E-2,1E2])
xlabel('Time step [s]')
ylabel('Position error norm after 3 hours[m]')
hold on
loglog(10.^integrationErrorBehaviour_long(:,1),sqrt(sum(integrationErrorBehaviour_long(:,2:4).^2')),'-*','LineWidth',2)

legend('64 bit (double) time and state scalar','80 bit (long double) state scalar and split time scalar')


set(gcf, 'Units', 'normalized', 'Position', [0,0,0.375 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 22.5 30]);
set(gcf,'PaperPositionMode','auto');

if( saveResults )
    saveas(gcf,strcat(saveFolder,'truncationRoundingErrorRk4'),'png');
end



