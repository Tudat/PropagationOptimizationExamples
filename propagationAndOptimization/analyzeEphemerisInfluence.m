clc
clear all
close all

currentPath = '/home/dominic/Software/tudatBundle/tudatApplications/propagationAndOptimization/';
dataPath = strcat(currentPath,'EphemerisSettings/');

orbit= cell(3,3)
energy= cell(3,3)

figure(1)
for i=1:3
    for j=1:3
        orbit{i,j}=load(strcat(dataPath,'stateMoonOrbiter_',num2str(i-1),'_',num2str(j-1),'.dat'));
        energy{i,j}=load(strcat(dataPath,'orbitalEnergyMoonOrbiter_',num2str(i-1),'_',num2str(j-1),'.dat'));
        if( i > 1 )
            subplot(3,3,(i-2)*3+j)
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{1,j}(:,2:4))
        end
        
        if( i == 3 )
            subplot(3,3,(i-1)*3+j)
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{2,j}(:,2:4))
        end       
    end
end
%%
figure(2)
semilogy(orbit{2,2}(:,1)/86400,orbit{3,2}(:,2:4)-orbit{1,2}(:,2:4),'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Position difference [m]')
legend('x','y','z','Location','NorthWest')
title('Orbiter state difference: DE - Kepler ephemerides')


