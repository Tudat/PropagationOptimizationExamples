clc
clear all
close all


saveFolder = '/home/dominic/Software/tudatBundle/tudat/Tudat/';

nominalState = load(strcat(saveFolder,'entryJ2SensitivityNominalTrajectory.dat'));
nominalPosition = sqrt(nominalState(:,2).^2+nominalState(:,3).^2+nominalState(:,4).^2);

nominalSize = max(size(nominalPosition));

numberOfSamples = 5000;
perturbedStates = cell(numberOfSamples,1);
statePerturbations = cell(numberOfSamples,1);
positionPerturbations = cell(numberOfSamples,1);
perturbedSizes = zeros(numberOfSamples,1);

for i=1:numberOfSamples
    i
    perturbedStates{i,1} = load(strcat(saveFolder,'entryJ2SensitivityPerturbedTrajectoryNoJ2',num2str(i-1),'.dat'));
    
    perturbedSizes(i) = max(size(perturbedStates{i,1}));
    numberOfRows = min([perturbedSizes(i) nominalSize]);
    
    statePerturbations{i,1} = perturbedStates{i,1}(1:numberOfRows,:) - nominalState(1:numberOfRows,:);
    
    positionPerturbations{i,1} =  sqrt(statePerturbations{i,1}(:,2).^2+statePerturbations{i,1}(:,3).^2+statePerturbations{i,1}(:,4).^2);
end
%%
close all
figure(1)
minimumSize = min(perturbedSizes);
for i=1:numberOfSamples
    plot(1:minimumSize,positionPerturbations{i,1}(1:minimumSize),'b')
    hold on
    finalError(i) = positionPerturbations{i,1}(minimumSize);

end
grid on
xlim([0 minimumSize])
xlabel('Time [s]')
ylabel('Position perturbation [m]')
%%
figure(2)
hist(finalError,25);

xlabel('Final position perturbation [m]')
ylabel('Number of occurences')