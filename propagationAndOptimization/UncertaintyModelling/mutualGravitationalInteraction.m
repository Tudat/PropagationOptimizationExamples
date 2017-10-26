clc
close all
clear all

preFitResidual = load(strcat(dataFolder,'moonMutualJ2PreFitResiduals.dat'));
postFitResidual = load(strcat(dataFolder,'moonMutualJ2PostFitResiduals.dat'));

numberOfEntries = max(size(preFitResidual));

preFitNorm = sqrt( preFitResidual(1:3:numberOfEntries).^2+preFitResidual(2:3:numberOfEntries).^2+preFitResidual(3:3:numberOfEntries).^2 );
postFitNorm = sqrt( postFitResidual(1:3:numberOfEntries).^2+postFitResidual(2:3:numberOfEntries).^2+postFitResidual(3:3:numberOfEntries).^2 );

timeVector = (1:numberOfEntries/3)*6*3600/86400;

figure(1)
subplot(1,2,1)
plot(timeVector,preFitResidual(1:3:numberOfEntries));
hold on
plot(timeVector,preFitResidual(2:3:numberOfEntries));
plot(timeVector,preFitResidual(3:3:numberOfEntries));
grid on
xlabel('Time [days]')
ylabel('Prefit position difference [m]')
legend('x','y','z','Location','NorthWest')

xlim([0 365])
subplot(1,2,2)
plot(timeVector,postFitResidual(1:3:numberOfEntries));
hold on
plot(timeVector,postFitResidual(2:3:numberOfEntries));
plot(timeVector,postFitResidual(3:3:numberOfEntries));
grid on
xlabel('Time [days]')
ylabel('Postfit position difference [m]')
xlim([0 365])

suptitle('Direct influence of lunar degree 2 gravity field on Moon position')

% plot(lunarOrbitPrefitResidual)
% hold on
% plot(lunarOrbitPostfitResidual)
