set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

set(0,'defaultAxesFontSize',10)
set(0,'defaultTextFontSize',10)

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
errorMap_Abam = cell(12,4,4);

for i=0
    for j=0:3
        functionEvaluations{i,j} = load(strcat(folder,'functionEvaluations_e_',num2str(i),'_intSett',num2str(j),fileSuffix,'.dat'));
        for k=9:13
            disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
            %errorMap{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),fileSuffix,'.dat'));
            errorMap_Abam{i+1,j+1,k-8}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),fileSuffix,'.dat'));
        end
    end
end

%%

set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontSize',12)

close all
eccentricities = [ 0.01, 0.05, 0.1, 0.25, 0.5, 0.9, 0.95, 0.99];
fixedStepSize = [2 4 8 16];
tolerances = [1E-15 1E-13 1E-11 1E-9];

titles = cell(5,1);
titles{1}='Variable step, variable order';
titles{2}='Variable step, order = 6';
titles{3}='Variable step, order = 8';
titles{4}='Variable step, order = 10';
titles{5}='Fixed step, variable order';


close all
for type = 0
    for i=1:7
        figure(i+8*type)
        for j=1:4
            for k=1:5
                subplot(2,3,k)
                maximumError(i,j,k) = max(sqrt(sum(errorMap_Abam{i,j,k}(:,2:4)'.^2)));
                semilogy(errorMap_Abam{i,j,k}(:,1),sqrt(sum(errorMap_Abam{i,j,k}(:,2:4)'.^2)),'LineWidth',2)
                
                hold on
                grid on
                
                if(j==4)
                    
                    xlim([0 14*86400])
                    
                    xlabel('t [s]');
                    ylabel('Position error [m]');
                    
                    title(titles{k})
                    if( j == 4 && k == 1 )
                        legend('tol=10E-15','tol=10E-13','tol=10E-11','tol=10E-9','Location','SouthEast')
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
    
end
%%
close all


for i=1:8
    figure(i+16)
    counter = 1;
    for j=1:4
        for k=1:5
            subplot(4,5,counter)
            sizes = size(errorMap_Abam{i,j,k}(:,1));
            numberOfTimeStep = sizes(1);
            semilogy(errorMap_Abam{i,j,k}(2:numberOfTimeStep,1),errorMap_Abam{i,j,k}(2:numberOfTimeStep,1)-errorMap_Abam{i,j,k}(1:(numberOfTimeStep-1),1));
            
            if( k ~= 5 )
                title(strcat(titles{k},' tol=',num2str(tolerances(j))))
            else
                title(strcat(titles{k},' dt=',num2str(fixedStepSize(j)),' s'))

            end
            hold on
            grid on
            
            counter = counter + 1;
            
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.5 0.5]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
    set(gcf,'PaperPositionMode','auto');
     if( saveResults )
         pause(1.0)
         saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_Abam',num2str(i)),'png');
     end
end
