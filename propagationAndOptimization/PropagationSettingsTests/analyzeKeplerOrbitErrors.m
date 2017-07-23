clc
close all
clear all

  set(0,'DefaultAxesFontSize',14)


errorMap = cell(12,4,4);
errorMap = cell(12,4,1);

for i=0:7
    for j=0:3
        for k=0:3
            disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k)))
            errorMap{i+1,j+1,k+1}=load(strcat('numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));
            if( k == 0 )
                errorBackwardsMap{i+1,j+1,k+1}=load(strcat('numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));

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
%%
close all
for i=1:8
    figure(i)
    for j=1:4
        for k=1:4
            subplot(2,2,k)
            if( ~(k==1 && j == 4 ) )
                semilogy(errorMap{i,j,k}(:,1),sqrt(sum(errorMap{i,j,k}(:,2:4)'.^2)),'LineWidth',2)
            end
            
            hold on
            grid on
            
            
            if(j==4)
                
                title(integratorTypes{k})
                
                ylim([1E-5 1E10])
                xlim([0 5*86400])
                
                xlabel('t [s]');
                ylabel('Position error [m]');
                
                if( k== 1)
                    legend('t=10 s','t=100 s', 't=1000 s')
                end
                
                if( k== 2)
                    legend('tol=10E-14','tol=10E-13','tol=10E-12','tol=10E-11')
                end
            end
        end
    end
    
    suptitle(strcat('Eccentricity=',num2str(eccentricities(i))));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
    saveas(figure(i),strcat('KeplerOrbitError',num2str(i)),'png');

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

%%
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
    saveas(figure(i+16),strcat('KeplerOrbitStepSizeControl',num2str(i)),'png');
end
%%
%%(i+1,j+1,k+1,l+1,m+1)
for i=0:2
    for j=0:3
        for k=0:3
            for l=0:1
                for m=1                    
                    disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(l)))
                    perturbedStateMap{i+1,j+1,k+1,l+1,2}=load(strcat('numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),'_accSett',num2str(m),'.dat'));
                    s=size(perturbedStateMap{i+1,j+1,k+1,l+1,2});
                    timeSteps{i+1,j+1,k+1,l+1,2}=perturbedStateMap{i+1,j+1,k+1,l+1,2}(2:s(1),1)-perturbedStateMap{i+1,j+1,k+1,l+1,2}((1:s(1)-1),1);
                    maximumTimeStep(i+1,j+1,k+1,l+1,m+1)=max(timeSteps{i+1,j+1,k+1,l+1,2});
                    minimumTimeStep(i+1,j+1,k+1,l+1,m+1)=min(timeSteps{i+1,j+1,k+1,l+1,2});
                    meanTimeStep(i+1,j+1,k+1,l+1,m+1)=mean(timeSteps{i+1,j+1,k+1,l+1,2});
                    standardDeviationTimeStep(i+1,j+1,k+1,l+1,m+1)=std(timeSteps{i+1,j+1,k+1,l+1,2});
                    numberOfTimeStepsPerturbed(i+1,j+1,k+1,l+1,m+1) = max(size(perturbedStateMap{i+1,j+1,k+1,l+1,2}));

                end
            end
        end
    end
end

for i=0:2
    for j=0:3
        for k=0:3
            for l=0
                for m=0                    
                    disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(l)))
                    perturbedStateMap{i+1,j+1,k+1,l+1,1}=load(strcat('numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),'_accSett',num2str(m),'.dat'));
                    s=size(perturbedStateMap{i+1,j+1,k+1,l+1,1});
                    timeSteps{i+1,j+1,k+1,l+1,1}=perturbedStateMap{i+1,j+1,k+1,l+1,1}(2:s(1),1)-perturbedStateMap{i+1,j+1,k+1,l+1,1}((1:s(1)-1),1);
                    maximumTimeStep(i+1,j+1,k+1,l+1,m+1)=max(timeSteps{i+1,j+1,k+1,l+1,1});
                    minimumTimeStep(i+1,j+1,k+1,l+1,m+1)=min(timeSteps{i+1,j+1,k+1,l+1,1});
                    meanTimeStep(i+1,j+1,k+1,l+1,m+1)=mean(timeSteps{i+1,j+1,k+1,l+1,1});
                    standardDeviationTimeStep(i+1,j+1,k+1,l+1,m+1)=std(timeSteps{i+1,j+1,k+1,l+1,1});
                    numberOfTimeStepsPerturbed(i+1,j+1,k+1,l+1,m+1) = max(size(perturbedStateMap{i+1,j+1,k+1,l+1,1}));
                end
            end
        end
    end
end

%%
close all

timeEvaluations = 0:10:3.0*3600;

i = 3;
k = 2;
m = 2;

figure(100)

for l=1:2
    referenceStateTimes = perturbedStateMap{i,1,k,l,m}(:,1);
    refereceStates = perturbedStateMap{i,1,k,l,m}(:,2:4);
    interpolatedRefernceStates = interp1(referenceStateTimes,refereceStates,timeEvaluations);
    subplot(1,2,l)
    for j=2:4
        currentStateTimes = perturbedStateMap{i,j,k,l,m}(:,1);
        currentStates = perturbedStateMap{i,j,k,l,m}(:,2:4);
        interpolatedStates = interp1(currentStateTimes,currentStates,timeEvaluations);

        scatter(timeEvaluations,sum(sqrt((interpolatedStates-interpolatedRefernceStates)'.^2)),'*');
        hold on
        grid on
        xlim([0 3*3600])
        
    end
    if( l == 1 )
       title('Cowell') 
    else
        title('Encke')
    end
    legend('Tol=10^{-13}','Tol=10^{-12}','Tol=10^{-11}');
    xlabel('Time [s]')
    ylabel('Position difference w.r.t. tol=10^{-14}')
        
end

 set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
saveas(figure(100),strcat('CowellEnckeComparePerturbed'),'png');



%%
load('integrationErrorBehaviour.dat');
close all

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
    saveas(gcf,'truncationRoundingErrorRk4','png');

grid on