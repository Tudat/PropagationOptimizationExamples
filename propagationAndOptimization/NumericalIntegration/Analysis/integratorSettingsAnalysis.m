set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

 close all
 clear all
 clc

folder = '../../SimulationOutput/NumericalIntegration/';
saveResults = false;

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
errorMap = cell(12,4,4);
numberOfEccentricities = 7;

for i=0:(numberOfEccentricities-1)
    for j=0:3
        tempEvaluations = load(strcat(folder,'functionEvaluations_e_',num2str(i),'_intSett',num2str(j),fileSuffix,'.dat'));
        functionEvaluations(i+1,j+1,:)= tempEvaluations(:,2)';
        
        for k=0:8
            disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
            
            errorMap{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),fileSuffix,'.dat'));
            if( k < 5 )
                errorMap_Abam{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k+9),'_intSett',num2str(j),fileSuffix,'.dat'));
            end
            if( k == 0 )
                errorBackwardsMap{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),fileSuffix,'.dat'));
                if( k < 6 )
                    %errorBackwardsMap_Abam{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k+9),'_intSett',num2str(j),fileSuffix,'.dat'));
                end
                
            end
        end
    end
end

%%
eccentricities = [ 0.01, 0.05, 0.1, 0.25, 0.5, 0.9, 0.95, 0.99];
fixedStepSize = [2 4 8 16];
abanmFixedStepSize = [2 4 8 16];

tolerances = toleranceMultiplier * [1E-15 1E-13 1E-11 1E-9];

integratorTypes = cell(4,1);
integratorTypes{1} = 'RK4';
integratorTypes{2} = 'RKF4(5)';
integratorTypes{3} = 'RKF5(6)';
integratorTypes{4} = 'RKF7(8)';

bsIntegratorTypes = cell(4,1);
bsIntegratorTypes{1} = 'BS4';
bsIntegratorTypes{2} = 'BS6';
bsIntegratorTypes{3} = 'BS8';
bsIntegratorTypes{4} = 'BS10';

abamTypes = cell(5,1);
abamTypes{1}='Variable step, variable order';
abamTypes{2}='Variable step, order = 6';
abamTypes{3}='Variable step, order = 8';
abamTypes{4}='Variable step, order = 10';
abamTypes{5}='Fixed step, variable order';
%%
close all
figureCounter = 1;
for type = 0:1
    for i=1:numberOfEccentricities
        figure(figureCounter)
        figureCounter = figureCounter + 1;
        for j=1:4
            for k=1:4
                subplot(2,2,k)
                if( ~(k==1 && j == 4 && type == 0 ) )
                    semilogy(errorMap{i,j,k+4*type}(:,1),sqrt(sum(errorMap{i,j,k+4*type}(:,2:4)'.^2)),'LineWidth',2)
                end
                if( type == 0 )
                    maximumError(i,j,k) = max(sqrt(sum(errorMap{i,j,k+4*type}(:,2:4)'.^2)));
                elseif( type == 1 )
                    maximumErrorBs(i,j,k) = max(sqrt(sum(errorMap{i,j,k+4*type}(:,2:4)'.^2)));
                end

                
                hold on
                grid on
                
                
                if(j==4)
                    
                    if( type == 0 )
                        title(strcat(integratorTypes{k},titleSuffix))
                    else
                        title(strcat(bsIntegratorTypes{k},titleSuffix))
                    end
                    
                    if( k== 1 && type == 0 )
                        ylim([1E-5 1E10])
                    else
                        ylim([1E-6 1E4])
                    end
                    
                    xlim([0 14*86400])
                    
                    xlabel('t [s]');
                    ylabel('Position error [m]');
                    
                    if( k== 1 && type == 0 )
                    legend('dt=2 s','dt=4 s','dt=8 s','dt=16 s')
                    end
                    
                    if( k== 2 || ( k == 1 && type == 1 ) )
                        if( useLongValues )
                            legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11')
                        else
                            legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9')
                        end
                    end
                end
            end
        end
        
        suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
        pause(0.1)
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
        set(gcf,'PaperPositionMode','auto');
        if( saveResults )
            pause(1.0)
            saveas(figure(i+8*type),strcat('KeplerOrbitError',num2str(i),'_',num2str(type),fileSuffix),'png');
        end
        
    end
    
end
%%
for i=1:numberOfEccentricities
    figure(figureCounter)
    figureCounter = figureCounter + 1;
    for j=1:4
        for k=1:5
            subplot(2,3,k)
            maximumErrorAbam(i,j,k) = max(sqrt(sum(errorMap_Abam{i,j,k}(:,2:4)'.^2)));
            semilogy(errorMap_Abam{i,j,k}(:,1),sqrt(sum(errorMap_Abam{i,j,k}(:,2:4)'.^2)),'LineWidth',2)
            
            hold on
            grid on
            
            if(j==4)
                
                xlim([0 14*86400])
                
                xlabel('t [s]');
                ylabel('Position error [m]');
                
                title(abamTypes{k})
                if( j == 4 && k == 1 )
                    if( useLongValues )
                        legend('tol=10E-17','tol=10E-15','tol=10E-13','tol=10E-11','Location','SouthEast')
                    else
                        legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','Location','SouthEast')
                    end
                elseif( j == 4 && k == 5 )
                    legend('dt=2 s','dt=4 s','dt=8 s','dt=16 s','Location','SouthEast')
                end
            end
        end
    end
    
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    pause(0.1)
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    
    if( saveResults )
        pause(1.0)
        saveas(figure(i+8*type),strcat('KeplerOrbitErrorAbam',num2str(i),'_',num2str(type),fileSuffix),'png');
    end
end
%%

clc
close all
lineStyles = cell(5);
lineStyles{1} = '-';
lineStyles{2} = '--';
lineStyles{3} = '-.';
lineStyles{4} = ':';


for i=1:numberOfEccentricities
    figure(figureCounter)
    for k=1:4
        loglog(squeeze(functionEvaluations(i,:,k)),squeeze(maximumError(i,1:4,k)),strcat('b',lineStyles{k},'*'))
                hold on
    end
    
    for k=1:4
        loglog(squeeze(functionEvaluations(i,:,k+4)),squeeze(maximumErrorBs(i,:,k)),strcat('r',lineStyles{k},'*'))
    end
    
    for k=1:4
        loglog(squeeze(functionEvaluations(i,:,k+8)),squeeze(maximumErrorAbam(i,:,k)),strcat('k',lineStyles{k},'*'))
    end
    grid on
    
    title(strcat('Eccentricity=',num2str(eccentricities(i))));
    legend('RK4','RK4(5)','RK5(6)','RK7(8)','BS4','BS6','BS8','BS10','ABAM','ABAM8','ABAM10','ABAM-fixed step','Location','NorthEastOutside')
    
    xlabel('Number of function evaluations [-]')
    ylabel('Maximum position error (2 week integration time) [m]')
    
%     for k=1:4
%         loglog(squeeze(functionEvaluations(i,1,k)),squeeze(maximumError(i,1,k)),strcat('b',lineStyles{k},'o'))
%     end
%     
%     for k=1:4
%         loglog(squeeze(functionEvaluations(i,1,k+4)),squeeze(maximumErrorBs(i,1,k)),strcat('r',lineStyles{k},'o'))
%     end
%     
%     for k=1:4
%         loglog(squeeze(functionEvaluations(i,1,k+8)),squeeze(maximumErrorAbam(i,1,k)),strcat('k',lineStyles{k},'o'))
%     end
    

    
    pause(0.1)
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
    set(gcf,'PaperPositionMode','auto');
    
    if( saveResults )
        saveas(gcf,strcat('accuracyTradeOff',num2str(i),'_',num2str(type),fileSuffix),'png');
    end
        figureCounter = figureCounter + 1;

end

grid on

%%
close all
for i=1:8
    for j=1:4
        figure(12+j)
        subplot(2,4,i)
        plot(errorMap{i,j,1}(:,1),sqrt(sum(errorMap{i,j,1}(:,2:4)'.^2)));
        hold on
        plot(errorBackwardsMap{i,j,1}(:,1),sqrt(sum(errorBackwardsMap{i,j,1}(:,2:4)'.^2)));
        grid on
        
        title(strcat('e=',num2str(eccentricities(i))));
        
        
        xlabel('t [s]');
        ylabel('Position error [m]');
        
        xlim([0 14*86400])
        
        if(i==8)
            legend('Forward','Backwards','Location','SouthEast')
        end
    end
    
end


for j=1:4
    figure(12+j)
    suptitle(strcat('Step size=',num2str(fixedStepSize(j)),'s'));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    
    if( saveResults )
        pause(1.0)
        saveas(figure(12+j),strcat('ForwardBackwardPropagation',num2str(j)),'png');
    end
    
end

%%
close all
for i=1:8
    figure(i+16)
    for j=1:4
        for k=2:4
            subplot(1,3,k-1)
            sizes = size(errorMap{i,j,k}(:,1));
            numberOfTimeStep = sizes(1)
            semilogy(errorMap{i,j,k}(2:numberOfTimeStep,1),errorMap{i,j,k}(2:numberOfTimeStep,1)-errorMap{i,j,k}(1:(numberOfTimeStep-1),1));
            
            hold on
            grid on
            
            if(j==4)
                
                title(integratorTypes{k})
                
                xlim([0 1.3E4])
                
                xlabel('t [s]');
                ylabel('Step size [s]');
                
                if( k== 1)
                    legend('t=10 s','t=100 s', 't=1000 s')
                end
                
                if( k== 2)
                    legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','Location','SouthEast')
                end
            end
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    if( saveResults )
        saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_RKF',num2str(i)),'png');
    end
end
%%
close all
for i=1:8
    figure(i+16)
    for j=1:4
        for k=2:4
            subplot(1,3,k-1)
            sizes = size(errorMap{i,j,k}(:,1));
            numberOfTimeStep = sizes(1)
            semilogy(errorMap{i,j,k}(2:numberOfTimeStep,1),errorMap{i,j,k}(2:numberOfTimeStep,1)-errorMap{i,j,k}(1:(numberOfTimeStep-1),1));
            
            hold on
            grid on
            
            if(j==4)
                
                title(integratorTypes{k})
                
                xlim([0 1.3E4])
                
                xlabel('t [s]');
                ylabel('Step size [s]');
                
                if( k== 1)
                    legend('t=10 s','t=100 s', 't=1000 s')
                end
                
                if( k== 2)
                    legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','Location','SouthEast')
                end
            end
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    if( saveResults )
        saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_RKF',num2str(i)),'png');
    end
end
%%
close all
for i=1:8
    figure(i+16)
    for j=1:4
        for k=1
            subplot(1,4,j)
            sizes = size(errorMap{i,j,k}(:,1));
            numberOfTimeStep = sizes(1)
            semilogy(errorMap{i,j,k}(2:numberOfTimeStep,1),errorMap{i,j,k}(2:numberOfTimeStep,1)-errorMap{i,j,k}(1:(numberOfTimeStep-1),1));
            
            hold on
            grid on
            
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
    set(gcf,'PaperPositionMode','auto');
    if( saveResults )
        pause(1.0)
        saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_BS',num2str(i)),'png');
    end
end
%%
i = 7;
k = 1;
j = 1;


keplerOrbitExample = load(strcat(folder,'numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett0_accSett0.dat'))
%%
close all
semilogy(keplerOrbitExample(:,1),sqrt(keplerOrbitExample(:,5).^2+keplerOrbitExample(:,6).^2+keplerOrbitExample(:,7).^2))
yyaxis right
semilogy(keplerOrbitExample(:,1),sqrt(keplerOrbitExample(:,2).^2+keplerOrbitExample(:,3).^2+keplerOrbitExample(:,4).^2))

xlim([0 15000])

grid on