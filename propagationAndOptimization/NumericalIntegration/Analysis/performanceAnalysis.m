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
numberOfEccentricities = 7;

for l=0
    for i=0:(numberOfEccentricities-1)
        for j=0:5
            tempEvaluations = load(strcat(folder,'functionEvaluations_e_',num2str(i),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'_accSett0.dat'));
            functionEvaluations_Rk(i+1,j+1,1:4,l+1)= tempEvaluations(1:4,2)';
            functionEvaluations_Bs(i+1,j+1,1:4,l+1)= tempEvaluations(5:8,2)';
            functionEvaluations_Rk(i+1,j+1,5,l+1)= tempEvaluations(9,2)';
            functionEvaluations_Abm(i+1,j+1,1:5,l+1)= tempEvaluations(10:14,2)';
            
            tempErrors =  load(strcat(folder,'forwardBackwardError_e_',num2str(i),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'_accSett0.dat'));
            forwardBackwardErrors_Rk(i+1,j+1,1:4,l+1)= tempErrors(1:4,3)';
            forwardBackwardErrors_Bs(i+1,j+1,1:4,l+1)= tempErrors(5:8,3)';
            forwardBackwardErrors_Rk(i+1,j+1,5,l+1)= tempErrors(9,3)';
            forwardBackwardErrors_Abm(i+1,j+1,1:5,l+1)= tempErrors(10:14,3)';
            
            for k=0:13
                disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
                if( k < 4 )
                    errorMap_Rk{i+1,j+1,k+1,l+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    errorMap_Rk_Back{i+1,j+1,k+1,l+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    maximumError_Rk(i+1,j+1,k+1,l+1) = max(errorMap_Rk{i+1,j+1,k+1,l+1}(:,2)');
                    
                elseif( k < 8 )
                    errorMap_Bs{i+1,j+1,k+1- 4,l+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    errorMap_Bs_Back{i+1,j+1,k+1-4,l+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    maximumError_Bs(i+1,j+1,k+1-4,l+1) = max(errorMap_Bs{i+1,j+1,k+1-4,l+1}(:,2)');
                elseif( k == 8 )
                    errorMap_Rk{i+1,j+1,5,l+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    errorMap_Rk_Back{i+1,j+1,5,l+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    maximumError_Rk(i+1,j+1,5,l+1) = max(errorMap_Rk{i+1,j+1,5,l+1}(:,2)');
                else
                    errorMap_Abm{i+1,j+1,k+1-9,l+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    errorMap_Abm_Back{i+1,j+1,k+1-9,l+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),fileSuffix,'.dat'));
                    maximumError_Abm(i+1,j+1,k+1-9,l+1) = max(errorMap_Rk{i+1,j+1,k+1-9,l+1}(:,2)');
                end
            end
        end
    end
end


%%

clc
close all

saveFolder = '/home/dominic/Documents/Courses/Propagation and Optimization/IntegrationFigures/';

eccentricities = [ 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95];
fixedStepSize = [2 4 8 16 32 64];

tolerances = toleranceMultiplier * [1E-15 1E-13 1E-11 1E-9];

propagatorTypes = cell(4,1);
propagatorTypes{1} = 'Cowell';
propagatorTypes{2} = 'Gauss-Kepler';
propagatorTypes{3} = 'Gauss-MEE';


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
for l = 1
    for addForwardBackwardLine = 0:1
        for i=1:numberOfEccentricities
            figure(figureCounter)
            figureCounter = figureCounter + 1;
            for j=1:6
                for k=1:5
                    subplot(2,3,k)
                    
                    if( k < 5 )
                        colorIndex = j;
                    else
                        colorIndex = 1 + j;
                    end
                    
                    semilogy(errorMap_Rk{i,j,k,l}(:,1),(errorMap_Rk{i,j,k,l}(:,2)'),'LineWidth',2,'Color',colors(colorIndex,:))
                    hold on
                    if( addForwardBackwardLine )
                        semilogy(errorMap_Rk{i,j,k,l}(:,1),ones(size(errorMap_Rk{i,j,k,l}(:,1)))*forwardBackwardErrors_Rk(i,j,k,l),'--','Color',colors(colorIndex,:))
                    end
                    
                    grid on
                    if(j==6 )
                        
                        if( k== 1)
                            legend('t=2 s','t=4 s','t=8 s','t=16 s','t=32 s','Location','SouthEast')
                        end
                        if( ~useLongValues )
                            
                            if( k== 2)
                                legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','tol=10E-5','Location','SouthEast')
                            end
                        else
                            
                            if( k== 2)
                                legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','Location','SouthEast')
                            end
                        end
                        
                        title(strcat(integratorTypes{k},titleSuffix))
                        
                        xlim([0 7*86400])
                        
                        xlabel('t [s]');
                        ylabel('Position error [m]');
                        
                    end
                end
            end
            
            suptitle(strcat('Eccentricity=',num2str(eccentricities(i)),{' '},propagatorTypes{l}));
            pause(0.1)
            set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
            set(gcf,'PaperPositionMode','auto');
            if( saveResults )
                pause(0.1)
                if( addForwardBackwardLine == 0 )
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_RK_',num2str(i),fileSuffix),'png');
                else
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_RK_withForwardBackward',num2str(i),fileSuffix),'png');
                end
            end
        end
        
        close all
        
        %%
        for i=1:numberOfEccentricities
            figure(figureCounter)
            figureCounter = figureCounter + 1;
            for j=1:6
                for k=1:4
                    subplot(2,2,k)
                    semilogy(errorMap_Bs{i,j,k,l}(:,1),(errorMap_Bs{i,j,k,l}(:,2)'),'LineWidth',2,'Color',colors(j,:))
                    hold on
                    if( addForwardBackwardLine )
                        semilogy(errorMap_Bs{i,j,k,l}(:,1),ones(size(errorMap_Bs{i,j,k,l}(:,1)))*forwardBackwardErrors_Bs(i,j,k,l),'--','Color',colors(j,:))
                    end
                    grid on
                    if(j==6)
                        title(strcat(bsIntegratorTypes{k},titleSuffix))
                        
                        xlim([0 7*86400])
                        
                        xlabel('t [s]');
                        ylabel('Position error [m]');
                        
                        if( ~useLongValues )
                            
                            if( k== 1)
                                legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','tol=10E-5','Location','SouthEast')
                            end
                        else
                            
                            if( k== 1)
                                legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','Location','SouthEast')
                            end
                        end
                        
                    end
                end
            end
            
            suptitle(strcat('Eccentricity=',num2str(eccentricities(i)),{' '},propagatorTypes{l}));
            pause(0.1)
            set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
            set(gcf,'PaperPositionMode','auto');
            
            if( saveResults )
                pause(0.1)
                if( addForwardBackwardLine == 0 )
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_Bs_',num2str(i),fileSuffix),'png');
                else
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_Bs_withForwardBackward',num2str(i),fileSuffix),'png');
                end
            end
        end
        
        %%
        close all
        
        for i=1:numberOfEccentricities
            figure(figureCounter)
            figureCounter = figureCounter + 1;
            for j=1:6
                for k=1:5
                    subplot(2,3,k)
                    semilogy(errorMap_Abm{i,j,k,l}(:,1),(errorMap_Abm{i,j,k,l}(:,2)'),'LineWidth',2,'Color',colors(j,:))
                    hold on
                    if( addForwardBackwardLine )
                        semilogy(errorMap_Abm{i,j,k,l}(:,1),ones(size(errorMap_Abm{i,j,k,l}(:,1)))*forwardBackwardErrors_Abm(i,j,k,l),'--','Color',colors(j,:))
                    end
                    grid on
                    if(j==6 )
                        title(strcat(abamTypes{k},titleSuffix))
                        
                        xlim([0 7*86400])
                        
                        xlabel('t [s]');
                        ylabel('Position error [m]');
                        
                        
                        if( ~useLongValues )
                            if( k== 1)
                                legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','tol=10E-5','Location','SouthEast')
                            end
                        else
                            if( k== 1)
                                legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','Location','SouthEast')
                            end
                        end
                        
                        if( k== 5)
                            legend('t=2 s','t=4 s','t=8 s','t=16 s','t=32 s','t=64 s')
                        end
                        
                    end
                end
            end
            
            suptitle(strcat('Eccentricity=',num2str(eccentricities(i)),{' '},propagatorTypes{l}));
            pause(0.1)
            set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
            set(gcf,'PaperPositionMode','auto');
            
            if( saveResults )
                pause(0.1)
                if( addForwardBackwardLine == 0 )
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_Abm_',num2str(i),fileSuffix),'png');
                else
                    saveas(figure(gcf),strcat(saveFolder,'KeplerOrbitError_Abm_withForwardBackward',num2str(i),fileSuffix),'png');
                end
            end
        end
    end
    
    close all
    
    
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
            loglog(squeeze(functionEvaluations_Rk(i,:,k,l)),squeeze(maximumError_Rk(i,:,k,l)),strcat('b',lineStyles{k},'*'),'LineWidth',lineWidths{k})
            hold on
        end
        
        for k=1:4
            loglog(squeeze(functionEvaluations_Bs(i,1:4,k,l)),squeeze(maximumError_Bs(i,1:4,k,l)),strcat('r',lineStyles{k},'*'),'LineWidth',lineWidths{k})
        end
        
        for k=1:5
            loglog(squeeze(functionEvaluations_Abm(i,1:5,k,l)),squeeze(maximumError_Abm(i,1:5,k,l)),strcat('k',lineStyles{k},'*'),'LineWidth',lineWidths{k})
        end
        grid on
        
        title(strcat('Eccentricity=',num2str(eccentricities(i))));
        legend('RK4','RK4(5)','RK5(6)','RK7(8)','DOPRI8(7)','BS4','BS6','BS8','BS10','ABM','ABM6','ABM8','ABM10','ABM-fixed step','Location','NorthEastOutside')
        
        xlabel('Number of function evaluations [-]')
        ylabel('Maximum position error (1 week integration time) [m]')
        
        pause(0.1)
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
        set(gcf,'PaperPositionMode','auto');
        
        if( saveResults )
            saveas(gcf,strcat(saveFolder,'accuracyTradeOff',num2str(i),fileSuffix),'png');
        end
        figureCounter = figureCounter + 1;
        
    end
    
    close all
    %%
    for i=1:numberOfEccentricities
        figure(figureCounter)
        figureCounter = figureCounter + 1;
        for j=1:6
            for k=2:5
                subplot(1,4,k-1)
                sizes = size(errorMap_Rk{i,j,k,l}(:,1));
                numberOfTimeStep = sizes(1);
                if( k < 5 )
                    colorIndex = j;
                else
                    colorIndex = 1 + j;
                end
                
                semilogy(errorMap_Rk{i,j,k,l}(2:numberOfTimeStep,1),errorMap_Rk{i,j,k,l}(2:numberOfTimeStep,1)-errorMap_Rk{i,j,k,l}(1:(numberOfTimeStep-1),1),'Color',colors(colorIndex,:));
                
                hold on
                grid on
                
                if(j==6)
                    title(integratorTypes{k})
                    xlim([0 1.3E4])
                    
                    xlabel('t [s]');
                    ylabel('Step size [s]');
                    
                    if( k== 1)
                        legend('t=2 s','t=4 s','t=8 s','t=16 s','t=32 s')
                    end
                    
                    if( ~useLongValues )
                        if( k== 2)
                            legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','tol=10E-5','Location','SouthEast')
                        end
                    else
                        if( k== 2)
                            legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','Location','SouthEast')
                        end
                    end
                end
            end
        end
        suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
        
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
        set(gcf,'PaperPositionMode','auto');
        if( saveResults )
            saveas(gcf,strcat(saveFolder,'KeplerOrbitStepSizeControl_RKF',num2str(i)),'png');
        end
    end
    %%
    close all
    
    for i=1:numberOfEccentricities
        figure(figureCounter)
        figureCounter = figureCounter + 1;
        for j=1:6
            for k=1:4
                subplot(1,4,k)
                sizes = size(errorMap_Bs{i,j,k,l}(:,1));
                numberOfTimeStep = sizes(1)
                semilogy(errorMap_Bs{i,j,k,l}(2:numberOfTimeStep,1),errorMap_Bs{i,j,k,l}(2:numberOfTimeStep,1)-errorMap_Bs{i,j,k,l}(1:(numberOfTimeStep-1),1),'Color',colors(j,:));
                
                hold on
                grid on
                
                if(j==6)
                    title(bsIntegratorTypes{k})
                    xlim([0 1.3E4])
                    
                    xlabel('t [s]');
                    ylabel('Step size [s]');
                    
                    if( ~useLongValues )
                        if( k== 1)
                            legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','tol=10E-5','Location','SouthEast')
                        end
                    else
                        if( k== 1)
                            legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','tol=10E-7','Location','SouthEast')
                        end
                    end
                    
                    if( k== 5)
                        legend('t=2 s','t=4 s','t=8 s','t=16 s','t=32 s','t=64 s')
                    end
                end
            end
        end
        suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
        
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
        set(gcf,'PaperPositionMode','auto');
        if( saveResults )
            saveas(gcf,strcat(saveFolder,'KeplerOrbitStepSizeControl_Bs',num2str(i)),'png');
        end
    end
end

