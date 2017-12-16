clc
close all
clear all


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
