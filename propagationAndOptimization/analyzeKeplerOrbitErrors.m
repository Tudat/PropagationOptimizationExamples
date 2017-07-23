clc
close all
clear all


errorMap = cell(12,4,4);
errorMap = cell(12,4,1);

for i=0:11
    for j=0:3
        for k=0:3
            errorMap{i+1,j+1,k+1}=load(strcat('numericalKeplerOrbitError_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));
            if( k == 0 )
                errorBackwardsMap{i+1,j+1,k+1}=load(strcat('numericalKeplerOrbitErrorBack_e_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'.dat'));

            end
        end
    end
end


for i=1:12
    figure(i)
    for j=1:4
        for k=1:4
            subplot(2,2,k)
            semilogy(errorMap{i,j,k}(:,1),sqrt(sum(errorMap{i,j,k}(:,2:4)'.^2)));
            hold on
            grid on
        end
    end
end
%%
for i=1:12
    for j=1:4
        figure(12+j)
        subplot(4,3,i)
        plot(errorMap{i,j,1}(:,1),sqrt(sum(errorMap{i,j,1}(:,2:4)'.^2)));
        hold on
        plot(errorBackwardsMap{i,j,1}(:,1),sqrt(sum(errorBackwardsMap{i,j,1}(:,2:4)'.^2)));
        grid on
        
    end
end

%%
close all
for i=1:12
    figure(i+16)
    for j=1:4
        for k=1:4
            subplot(2,2,k)
            sizes = size(errorMap{i,j,k}(:,1));
            numberOfTimeStep = sizes(1)
            plot(errorMap{i,j,k}(2:numberOfTimeStep,1),errorMap{i,j,k}(2:numberOfTimeStep,1)-errorMap{i,j,k}(1:(numberOfTimeStep-1),1));
            hold on
            grid on
        end
    end
end

%%
perturbedStateMap=cell(3,3,3,2)
for i=0:2
    for j=1:3
        for k=1:3
            for l=0:1
                perturbedStateMap{i+1,j,k,l+1}=load(strcat('numericalKeplerOrbit_eccSett_',num2str(i),'_intType',num2str(k),'_intSett',num2str(j),'_propSett',num2str(k),'_accSett1.dat'));
            end
        end
    end
end