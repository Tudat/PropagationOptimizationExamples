clc
clear all
close all

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

dataPath = '../../SimulationOutput/AccelerationModels/';

orbit= cell(3,3);
orbit_full= cell(3,3);

figure(1)
for i=1:3
    for j=1:3
        orbit{i,j}=load(strcat(dataPath,'stateMoonOrbiterOriginCases_',num2str(i-1),'_',num2str(j-1),'.dat'));
        
        if( j > 1 )
            subplot(3,3,(i-1)*3+(j-1))
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{i,1}(:,2:4))
        end
        
        if( j == 3 )
          subplot(3,3,(i-1)*3+(j))
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{i,2}(:,2:4))  
        end
    end
end

figure(2)
for i=1:3
    for j=1:3
        orbit{i,j}=load(strcat(dataPath,'stateMoonOrbiterOriginCases_',num2str(i-1),'_',num2str(j-1),'.dat'));
        
        if( i > 1 )
            subplot(3,3,(i-2)*3+(j))
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{1,j}(:,2:4))
        end
        
        if( i == 3 )
             subplot(3,3,(i-1)*3+(j))
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{2,j}(:,2:4))
        end
    end
end

%%
close all
figure(3)
subplot(1,2,1)
plot(orbit{2,2}(:,1)/86400,orbit{1,2}(:,2:4)-orbit{1,1}(:,2:4),'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Position difference [m]')
legend('x','y','z','Location','NorthWest')
title('Earth and Sun')
subplot(1,2,2)
plot(orbit{2,2}(:,1)/86400,orbit{1,3}(:,2:4)-orbit{1,2}(:,2:4),'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Position difference [m]')
title('Jup., Mars and Venus')

suptitle('Effect of third-body perturbations, Moon-centered propagation')

figure(4)
subplot(1,2,1)
plot(orbit{2,2}(:,1)/86400,orbit{3,2}(:,2:4)-orbit{3,1}(:,2:4),'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Position difference [m]')
legend('x','y','z','Location','NorthWest')
title('Earth and Sun')

subplot(1,2,2)
plot(orbit{2,2}(:,1)/86400,orbit{3,3}(:,2:4)-orbit{3,2}(:,2:4),'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Position difference [m]')
legend('x','y','z','Location','NorthWest')
title('Jup., Mars and Venus')
suptitle('Effect of third-body perturbations, Sun-centered propagation')

figure(5)
subplot(1,2,1)
plot3(orbit{1,1}(:,2),orbit{1,1}(:,3),orbit{1,1}(:,4))
grid on
xlabel('x-pos. [m]')
ylabel('y-pos. [m]')
zlabel('z-pos. [m]')
title('Moon-centered')

subplot(1,2,2)
plot3(orbit{3,1}(:,2),orbit{3,1}(:,3),orbit{3,1}(:,4))
grid on
xlabel('x-pos. [m]')
ylabel('y-pos. [m]')
zlabel('z-pos. [m]')
title('Sun-centered')

suptitle('Moon-centered orbit, using only Moon point mass acceleration')



for figureIndex = 3:5
    set(figure(figureIndex), 'Units', 'normalized', 'Position', [0,0,0.5,0.5]);
    set(figure(figureIndex),'PaperUnits','centimeters','PaperPosition',[0 0 30 15 ])
    set(figure(figureIndex),'PaperPositionMode','auto')
end
%%
saveas(figure(3),strcat('PerturbationInfluence_MoonCentered'),'png')
saveas(figure(4),strcat('PerturbationInfluence_SunCentered'),'png')
saveas(figure(5),strcat('MoonOrbits_MoonAndSunCentered'),'png')

%%
sunState=load(strcat(dataPath,'sunStateComparison.dat'));
earthState=load(strcat(dataPath,'earthStateComparison.dat'));
moonState=load(strcat(dataPath,'moonStateComparison.dat'));
%%
close all
correctOrbit = orbit{1,1}
wrongOrbit = orbit{2,1};

figure(6)
subplot(1,2,1)
plot(correctOrbit(:,1),correctOrbit(:,2),'b')
hold on
plot(correctOrbit(:,1),correctOrbit(:,3),'r')

subplot(1,2,2)
plot(correctOrbit(:,1),wrongOrbit(:,2),'b')
hold on
plot(correctOrbit(:,1),wrongOrbit(:,3),'r')
