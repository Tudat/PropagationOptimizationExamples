set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontSize',12)

clc
close all
clear all

saveFolder = '/home/dominic/Documents/Courses/Propagation and Optimization/IntegrationFigures/';


for j=1:4
    
    startValue = 1;
    step=10^(-10+j);
    trueValue = startValue + step;
    totalNumberOfSteps = 100;
    finalValueError = zeros(totalNumberOfSteps,1);
    
    for test = 1:totalNumberOfSteps
        currentValue = startValue;
        numberOfSteps = 1E3 * test;
        extraTerm = step / (1E3 * test);
        
        for i=1:numberOfSteps
            currentValue = currentValue + extraTerm;
        end
        
        finalValue(test,1) = currentValue - trueValue;
        finalValueError(test,1) = currentValue - trueValue;
        numberOfStepsList(test,1) = numberOfSteps;
    end
    
    trueValues(j,1) = trueValue;
    
    
    subplot(2,2,j)
    
    scatter(numberOfStepsList,finalValueError,'*b')
    grid on
    xlabel('N [-]','interpreter','latex')
    ylabel('Error [-]','interpreter','latex')
    
    title(strcat('Error in computation of $1 + N\alpha$ ($N\alpha=$',num2str(step),')'),'interpreter','latex')
    
end

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.5]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

saveas(gcf,strcat(saveFolder,'roundingErrorExample'),'png');
