set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',20)
set(0,'defaultTextFontSize',20)

close all
clear all
clc

dataFolder = '../../SimulationOutput/EquationsOfMotion/';

perturbedSatellite = load(strcat(dataFolder,'singlePerturbedSatellitePropagationHistory.dat'));
unperturbedSatellite = load(strcat(dataFolder,'singleUnperturbedSatellitePropagationKeplerianHistory.dat'));
s = max(size(unperturbedSatellite));

time = ( perturbedSatellite(:,1) - perturbedSatellite(1,1) )/3600;
derivativeTime = ( perturbedSatellite(2:s,1) - perturbedSatellite(2,1) )/3600;

unperturbedSatelliteDerivative = ( unperturbedSatellite(2:s,:) - unperturbedSatellite(1:(s-1),:) )/10;
perturbedSatelliteDerivative = ( perturbedSatellite(2:s,:) - perturbedSatellite(1:(s-1),:) )/10;
%%
close all
figure(1)

subplot(1,2,1)
plot(time,perturbedSatellite(:,2),'LineWidth',2)
hold on
plot(time,perturbedSatellite(:,3),'LineWidth',2)
plot(time,perturbedSatellite(:,4),'LineWidth',2)
grid on
legend('x','y','z')
xlabel('t [hr]')
ylabel('Position [m]')
xlim([0 24])

subplot(1,2,2)
xlim([0 24])

subplot(1,2,2)
plot(time,perturbedSatellite(:,5),'LineWidth',2)
hold on
plot(time,perturbedSatellite(:,6),'LineWidth',2)
plot(time,perturbedSatellite(:,7),'LineWidth',2)
grid on
xlabel('t [hr]')
ylabel('Velocity [m/s]')
xlim([0 24])

suptitle('Earth orbiter - Cowell Propagator')

set( figure(1), 'Units', 'normalized', 'Position', [0,0,1,0.75]);
set( figure(1),'PaperUnits','centimeters','PaperPosition',[0 0 60 30]);

pause(2.0)

saveas(gcf,strcat('exampleCartesianPropagation'),'epsc');
%%
figure(10)

plot(derivativeTime,perturbedSatelliteDerivative(:,5),'LineWidth',2)

hold on
plot(derivativeTime,perturbedSatelliteDerivative(:,6),'LineWidth',2)
plot(derivativeTime,perturbedSatelliteDerivative(:,7),'LineWidth',2)
grid on
xlabel('t [hr]')
ylabel('Total Acceleration [m/s$^2$]')
legend('x','y','z')
xlim([0 24])

set( figure(10), 'Units', 'normalized', 'Position', [0,0,0.5,0.75]);
set( figure(10),'PaperUnits','centimeters','PaperPosition',[0 0 30 30]);

pause(2.0)

saveas(gcf,strcat('cowellDerivative'),'epsc');
figure(11)

plot(derivativeTime,perturbedSatelliteDerivative(:,5)-unperturbedSatelliteDerivative(:,5),'LineWidth',2)

hold on
plot(derivativeTime,perturbedSatelliteDerivative(:,6)-unperturbedSatelliteDerivative(:,6),'LineWidth',2)
plot(derivativeTime,perturbedSatelliteDerivative(:,7)-unperturbedSatelliteDerivative(:,7),'LineWidth',2)
grid on
xlabel('t [hr]')
ylabel('Perturbing Acceleration [m/s$^2$]')
legend('x','y','z')
xlim([0 24])


set( figure(11), 'Units', 'normalized', 'Position', [0,0,0.5,0.75]);
set( figure(11),'PaperUnits','centimeters','PaperPosition',[0 0 30 30]);
pause(2.0)

saveas(gcf,strcat('enckeDerivative'),'epsc');
%%

figure(2)

subplot(1,2,1)
plot(time,perturbedSatellite(:,2)-unperturbedSatellite(:,2),'LineWidth',2)
hold on
plot(time,perturbedSatellite(:,3)-unperturbedSatellite(:,3),'LineWidth',2)
plot(time,perturbedSatellite(:,4)-unperturbedSatellite(:,4),'LineWidth',2)
grid on
legend('x','y','z')
xlabel('t [hr]')
ylabel('Position [m]')
xlim([0 24])

subplot(1,2,2)
plot(time,perturbedSatellite(:,5)-unperturbedSatellite(:,5),'LineWidth',2)
hold on
plot(time,perturbedSatellite(:,6)-unperturbedSatellite(:,6),'LineWidth',2)
plot(time,perturbedSatellite(:,7)-unperturbedSatellite(:,7),'LineWidth',2)
grid on
xlabel('t [hr]')
ylabel('Velocity [m/s]')
xlim([0 24])

suptitle('Earth orbiter - Encke Propagator')

set( figure(2), 'Units', 'normalized', 'Position', [0,0,1,0.75]);
set( figure(2),'PaperUnits','centimeters','PaperPosition',[0 0 60 30]);

pause(2.0)

saveas(gcf,strcat('exampleCartesianPropagationDifference'),'epsc');
%%
perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationKeplerianHistory.dat'));

figure(3)
for i=1:6
    if( i > 2 )
        multiplier = 180/pi;
    else
        multiplier = 1;
    end
    subplot(2,3,i)
    plot(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
    xlim([0 24])
    xlabel('t [hr]')
    grid on
    if( i == 1)
        ylabel('a [m]')
    elseif( i == 2)
        ylabel('e [-]')
    elseif( i == 3)
        ylabel('i [deg]')
    elseif( i == 4)
        ylabel('$\omega$ [deg]')
    elseif( i == 5)
        ylabel('$\Omega$ [deg]')
    elseif( i == 6)
        ylabel('$\thet$a [deg]')
    end
end

suptitle('Earth orbiter - Kepler Propagator')

set( figure(3), 'Units', 'normalized', 'Position', [0,0,1,0.75]);
set( figure(3),'PaperUnits','centimeters','PaperPosition',[0 0 60 30]);

pause(2.0)

saveas(gcf,strcat('exampleKeplerPropagation'),'epsc');
%%
close all
perturbedSatelliteKeplerDerivative = ( perturbedSatelliteKepler(2:s,:) - perturbedSatelliteKepler(1:(s-1),:) )/10;
perturbedSatelliteKeplerDerivative(:,7) = perturbedSatelliteKeplerDerivative(:,7) +(perturbedSatelliteKeplerDerivative(:,7)<0).*2*pi/10;
figure(5)
for i=1:6
    if( i > 2 )
        multiplier = 180/pi;
    else
        multiplier = 1;
    end
    subplot(2,3,i)
    plot(derivativeTime,multiplier*perturbedSatelliteKeplerDerivative(:,i+1),'LineWidth',2)
    xlim([0 24])
    xlabel('t [hr]')
    grid on
    if( i == 1)
        ylabel('$da/dt$ [m/s]','interpreter','latex')
    elseif( i == 2)
        ylabel('$de/dt$ [1/s]','interpreter','latex')
    elseif( i == 3)
        ylabel('$di/dt$ [deg/s]','interpreter','latex')
    elseif( i == 4)
        ylabel('$d\omega/dt$ [deg/s]','interpreter','latex')
    elseif( i == 5)
        ylabel('$d\Omega/dt$ [deg/s]','interpreter','latex')
    elseif( i == 6)
        ylabel('$d\theta/dt$ [deg/s]','interpreter','latex')
    end
end

suptitle('Earth orbiter - Kepler Element Derivatives')

set( figure(5), 'Units', 'normalized', 'Position', [0,0,1,0.75]);
set( figure(5),'PaperUnits','centimeters','PaperPosition',[0 0 60 30]);

pause(2.0)

saveas(gcf,strcat('exampleKeplerDerivatives'),'epsc');

%%
perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationMeeHistory.dat'));

figure(4)
for i=1:6
        if( i == 6 )
            multiplier = 180/pi;
        else
            multiplier = 1;
        end
    multiplier = 1;
    subplot(2,3,i)
    plot(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
    xlim([0 24])
    xlabel('t [hr]')
    grid on
    if( i == 1)
        ylabel('p [m]')
    elseif( i == 2)
        ylabel('f [-]')
    elseif( i == 3)
        ylabel('g [-]')
    elseif( i == 4)
        ylabel('h [-]')
    elseif( i == 5)
        ylabel('k [-]')
    elseif( i == 6)
        ylabel('L [deg]')
    end
end

suptitle('Earth orbiter - MEE Propagator')

set( figure(4), 'Units', 'normalized', 'Position', [0,0,1,0.75]);
set( figure(4),'PaperUnits','centimeters','PaperPosition',[0 0 60 30]);

pause(2.0)

saveas(gcf,strcat('exampleMeePropagation'),'epsc');

%%

rRatio = 0.75:0.001:1.25;
q=0.5*(rRatio.^2-1);
qFunc=2*q./(1+2*q).*(1+1./(1+2*q+sqrt(1+2*q)));
close all
plot(rRatio,qFunc,'LineWidth',2)
grid on
ylabel('$\mathcal{F}(q)$ [-]','interpreter','latex')
xlabel('$r/ \rho [-]$','interpreter','latex')
xlim([0.75 1.25])

title('q-function for Encke propagator')

set( gcf, 'Units', 'normalized', 'Position', [0,0,0.5,0.75]);
set( gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 30]);
pause(2.0)
saveas(gcf,strcat('qFunction'),'epsc');
