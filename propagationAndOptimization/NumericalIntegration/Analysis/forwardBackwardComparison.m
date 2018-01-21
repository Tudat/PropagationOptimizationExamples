set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

close all
clear all
clc

folder = '../../SimulationOutput/NumericalIntegration/';
saveResults = true;

useLongValues = true;

if( useLongValues )
    fileSuffix = '_long' ;
    titleSuffix = ', long double scalars';
    toleranceMultiplier = 0.01;
else
    fileSuffix = '' ;
    titleSuffix = '';
    toleranceMultiplier = 1.0;
end
errorMap = cell(12,4,4);
numberOfEccentricities = 8;

for i=0
    for j=0:1                
        for k=0:1
           
            forwardPropagation{i+1,j+1,k+1}=load(strcat(folder,'perturbedOrbit_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(k),'_propSett1.dat'));
            backwardPropagation{i+1,j+1,k+1}=load(strcat(folder,'perturbedOrbitBackward_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(k),'_propSett1.dat'));

        end
    end
end

%%
clc
plot(forwardPropagation{1,1,1}(:,1),forwardPropagation{1,1,1}(:,2))
hold on
plot(backwardPropagation{1,1,1}(:,1),backwardPropagation{1,1,1}(:,2),'r--')
