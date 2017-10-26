clc
close all
clear all

  set(0,'DefaultAxesFontSize',14)


errorMap = cell(12,4,4);
errorMap = cell(12,4,1);

folder = '../../SimulationOutput/IntergratorAndPropagatorInfluence/';
for i=0:7
    for j=0:3
        for k=0:7
            disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
            errorMap{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));
            if( k == 0 )
                errorBackwardsMap{i+1,j+1,k+1}=load(strcat(folder,'numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));

            end
        end
    end
end

%%
eccentricities = [ 0.01, 0.05, 0.1, 0.25, 0.5, 0.9, 0.95, 0.99];
fixedStepSize = [10 100 1000 10000];
tolerances = [1E-14 1E-13 1E-12 1E-11];

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
%%
close all
for type = 0:1
    for i=1:8
        figure(i+8*type)
        for j=1:4
            for k=1:4
                subplot(2,2,k)
                if( ~(k==1 && j == 4 && type == 0 ) )
                    semilogy(errorMap{i,j,k+4*type}(:,1),sqrt(sum(errorMap{i,j,k+4*type}(:,2:4)'.^2)),'LineWidth',2)
                end
                
                hold on
                grid on
                
                
                if(j==4)
                    
                    if( type == 0 )
                        title(integratorTypes{k})
                    else
                        title(bsIntegratorTypes{k})
                    end
                    
                    if( type == 0 )
                        ylim([1E-5 1E10])
                    else
                        ylim([1E-6 1E4])
                    end
                    
                    xlim([0 14*86400])
                    
                    xlabel('t [s]');
                    ylabel('Position error [m]');
                    
                    if( k== 1 && type == 0 )
                        legend('t=10 s','t=100 s', 't=1000 s')
                    end
                    
                    if( k== 2 || ( k == 1 && type == 1 ) )
                        legend('tol=10E-14','tol=10E-13','tol=10E-12','tol=10E-11')
                    end
                end
            end
        end
        
        suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
        pause(0.1)
        set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
        set(gcf,'PaperPositionMode','auto');
        pause(1.0)
        saveas(figure(i),strcat('KeplerOrbitError',num2str(i),'_',num2str(type)),'png');
        
    end
    
end
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
        
        xlim([0 5*86400])
        
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
    saveas(figure(12+j),strcat('ForwardBackwardPropagation',num2str(j)),'png');

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
                    legend('tol=10E-14','tol=10E-13','tol=10E-12','tol=10E-11','Location','SouthEast')
                end
            end
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_RKF',num2str(i)),'png');
end

%%
close all
for i=1:8
    figure(i+16)
    for j=1:4
        for k=1:4
            subplot(1,4,k)
            sizes = size(errorMap{i,j,k+4}(:,1));
            numberOfTimeStep = sizes(1)
            semilogy(errorMap{i,j,k+4}(2:numberOfTimeStep,1),errorMap{i,j,k+4}(2:numberOfTimeStep,1)-errorMap{i,j,k+4}(1:(numberOfTimeStep-1),1));
            
            hold on
            grid on
            
            if(j==4)
                
                title(bsIntegratorTypes{k})
                
                xlim([0 1.3E4])
                ylim([1.0E-3 1.0E4])

                xlabel('t [s]');
                ylabel('Step size [s]');

                
                if( k== 1)
                    legend('tol=10E-14','tol=10E-13','tol=10E-12','tol=10E-11','Location','SouthEast')
                end
            end
        end
    end
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl_BS',num2str(i)),'png');
end