set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

close all
clear all
clc

folder = '../../SimulationOutput/IntergratorAndPropagatorInfluence/';
saveResults = true;

for i=0:2
    for j=0:3
        for k=0:3
            for l=0:1
                for m=1
                    disp(strcat(num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(l)))
                    perturbedStateMap{i+1,j+1,k+1,l+1,2}=load(strcat(folder,'numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),'_accSett',num2str(m),'.dat'));
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
                    perturbedStateMap{i+1,j+1,k+1,l+1,1}=load(strcat(folder,'numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(l),'_accSett',num2str(m),'.dat'));
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

if( saveResults )
    saveas(figure(100),strcat('CowellEnckeComparePerturbed'),'png');
end