clc
clear all
close all

saveFolder = '/home/dominic/Software/tudatBundle/tudat/bin/unit_tests/';
preFitResidual = load(strcat(saveFolder,'ApolloPreFitResiduals.dat'));
postFitResidual = load(strcat(saveFolder,'ApolloPostFitResiduals.dat'));

close all
numberOfEntries = max(size(postFitResidual));

preFitNorm = sqrt( preFitResidual(1:3:numberOfEntries).^2+preFitResidual(2:3:numberOfEntries).^2+preFitResidual(3:3:numberOfEntries).^2 );
postFitNorm = sqrt( postFitResidual(1:3:numberOfEntries).^2+postFitResidual(2:3:numberOfEntries).^2+postFitResidual(3:3:numberOfEntries).^2 );

figure(1)
semilogy(1:numberOfEntries/3,preFitNorm);
hold on
semilogy(1:numberOfEntries/3,postFitNorm);

legend('Prefit','Postfit','Location','NorthWest')
grid on

time = (1:numberOfEntries/3);

figure(2)
subplot(1,2,1)
for i=1:3
scatter(time,preFitResidual(i:3:numberOfEntries),'.');
hold on
end
grid on

xlabel('Time [s]')
ylabel('Position difference [m]')
legend('x','y','z','Location','NorthWest')
title('Pre-fit influence')


subplot(1,2,2)
for i=1:3
scatter(time,postFitResidual(i:3:numberOfEntries),'.');
hold on
end
grid on

xlabel('Time [s]')
ylabel('Position difference [m]')
title('Post-fit influence')





