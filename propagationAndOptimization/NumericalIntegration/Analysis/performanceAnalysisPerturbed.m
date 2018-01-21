set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

close all
clear all
clc

folder = '../../SimulationOutput/NumericalIntegration/';
saveResults = true;

useLongValues = false;

if( useLongValues )
    fileSuffix = '_long' ;
    titleSuffix = ', long double scalars';
    toleranceMultiplier = 0.01;
else
    fileSuffix = '' ;
    titleSuffix = '';
    toleranceMultiplier = 1.0;
end
numberOfEccentricities = 3;
 
% Stopped at Forward/backward 0 0 604800   -0.113191  -0.0602384   -0.119346 0.000264293 0.000121275 8.64794e-05 0 0 0 0 0 0
% 3 0 1 1
for l=0:3
    for i=0:(numberOfEccentricities-1)
        for j=0:5
            tempEvaluations = load(strcat(folder,'functionEvaluations_e_',num2str(i),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'_accSett1.dat'));
            functionEvaluations_Rk(i+1,j+1,1:4,l+1)= tempEvaluations(1:4,2)';
            functionEvaluations_Bs(i+1,j+1,1:4,l+1)= tempEvaluations(5:8,2)';
            functionEvaluations_Rk(i+1,j+1,5,l+1)= tempEvaluations(9,2)';
            functionEvaluations_Abm(i+1,j+1,1:5,l+1)= tempEvaluations(10:14,2)';
            
            tempErrors =  load(strcat(folder,'forwardBackwardError_e_',num2str(i),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'_accSett1.dat'));
            forwardBackwardErrors_Rk(i+1,j+1,1:4,l+1)= tempErrors(1:4,3)';
            forwardBackwardErrors_Bs(i+1,j+1,1:4,l+1)= tempErrors(5:8,3)';
            forwardBackwardErrors_Rk(i+1,j+1,5,l+1)= tempErrors(9,3)';
            forwardBackwardErrors_Abm(i+1,j+1,1:5,l+1)= tempErrors(10:14,3)';
%             for k=0:14
%                 disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
%                 if( j > 0 )
%                     if( k < 4 )
%                         errorMap_Rk{i+1,j,k+1,l+1}=load(strcat(folder,'interpolatedPerturbedOrbit_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                         errorMap_Rk_B{i+1,j,k+1,l+1}=load(strcat(folder,'interpolatedPerturbedOrbitB_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                     elseif( k < 8 )
%                         errorMap_Bs{i+1,j,k+1- 4,l+1}=load(strcat(folder,'interpolatedPerturbedOrbit_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                         errorMap_Bs_B{i+1,j,k+1- 4,l+1}=load(strcat(folder,'interpolatedPerturbedOrbitB_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                         
%                     elseif( k == 8 )
%                         errorMap_Rk{i+1,j,5,l+1}=load(strcat(folder,'interpolatedPerturbedOrbit_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                         errorMap_Rk_B{i+1,j,5,l+1}=load(strcat(folder,'interpolatedPerturbedOrbitB_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                     else
%                         errorMap_Abm{i+1,j,k+1-9,l+1}=load(strcat(folder,'interpolatedPerturbedOrbit_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                         errorMap_Abm_B{i+1,j,k+1-9,l+1}=load(strcat(folder,'interpolatedPerturbedOrbitB_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
%                     end
%                 end
%             end
        end
        
    end
end

save('processedResults_Perturbed')



%%

clc
close all

saveFolder = '/home/dominic/Documents/Courses/Propagation and Optimization/IntegrationFigures/';

eccentricities = [ 0.01, 0.1, 0.5, 0.9];
fixedStepSize = [2 4 8 16 32 64];

tolerances = toleranceMultiplier * [1E-15 1E-13 1E-11 1E-9];

propagatorTypes = cell(4,1);
propagatorTypes{1} = 'Cowell';
propagatorTypes{2} = 'Gauss-Kepler';
propagatorTypes{3} = 'Gauss-MEE';
propagatorTypes{4} = 'Encke';


integratorTypes = cell(4,1);
integratorTypes{1} = 'RK4';
integratorTypes{2} = 'RKF4(5)';
integratorTypes{3} = 'RKF5(6)';
integratorTypes{4} = 'RKF7(8)';
integratorTypes{5} = 'DOPRI8(7)';

bsIntegratorTypes = cell(4,1);
bsIntegratorTypes{1} = 'BS4';
bsIntegratorTypes{2} = 'BS6';
bsIntegratorTypes{3} = 'BS8';
bsIntegratorTypes{4} = 'BS10';

abamTypes = cell(5,1);
abamTypes{1}='ABM var. step, var. order';
abamTypes{2}='ABM var. step, order = 6';
abamTypes{3}='ABM var. step, order = 8';
abamTypes{4}='ABM var. step, order = 10';
abamTypes{5}='ABM fixed step, var. order';

set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontSize',12)

colors = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

close all
figureCounter = 1;
l = 1;
%%
close all
for l = 1:4
    
    
    
    lineStyles = cell(5);
    lineStyles{1} = '-';
    lineStyles{2} = '--';
    lineStyles{3} = '-.';
    lineStyles{4} = ':';
    lineStyles{5} = ':';
    
    lineWidths = cell(5);
    lineWidths{1} = 1;
    lineWidths{2} = 1;
    lineWidths{3} = 1;
    lineWidths{4} = 1;
    lineWidths{5} = 2;
    
    for i=1:numberOfEccentricities
        figure(figureCounter)
        for k=1:5
            loglog(squeeze(functionEvaluations_Rk(i,:,k,l)),squeeze(forwardBackwardErrors_Rk(i,:,k,l)),strcat('b',lineStyles{k},'*'),'LineWidth',lineWidths{k})
            hold on
        end
        
        for k=1:4
            loglog(squeeze(functionEvaluations_Bs(i,1:4,k,l)),squeeze(forwardBackwardErrors_Bs(i,1:4,k,l)),strcat('r',lineStyles{k},'*'),'LineWidth',lineWidths{k})
        end
        
        for k=1:5
            loglog(squeeze(functionEvaluations_Abm(i,1:5,k,l)),squeeze(forwardBackwardErrors_Abm(i,1:5,k,l)),strcat('k',lineStyles{k},'*'),'LineWidth',lineWidths{k})
        end
        grid on
        
        title(strcat('Eccentricity=',num2str(eccentricities(i)),{' '},propagatorTypes{l}));
        legend('RK4','RK4(5)','RK5(6)','RK7(8)','DOPRI8(7)','BS4','BS6','BS8','BS10','ABM','ABM6','ABM8','ABM10','ABM-fixed step','Location','NorthEastOutside')
        
        xlabel('Number of function evaluations [-]')
        ylabel('Maximum position error (1 week integration time) [m]')
        
        ylim([1.0E-6 1.0E8])
        xlim([1E2 1E7])
        
        pause(0.1)
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
        set(gcf,'PaperPositionMode','auto');
        
        if( saveResults )
            saveas(gcf,strcat(saveFolder,'accuracyTradeOffPerturbed',num2str(i),'_',num2str(l),fileSuffix),'png');
        end
        figureCounter = figureCounter + 1;
        
    end
    
end

