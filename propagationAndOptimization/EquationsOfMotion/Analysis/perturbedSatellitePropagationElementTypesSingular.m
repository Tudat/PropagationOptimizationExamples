set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

set(0,'defaultAxesFontSize',20)
set(0,'defaultTextFontSize',20)

close all
clear all
clc

for k=1:3
    dataFolder = '../../SimulationOutput/EquationsOfMotion/';
    perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationKeplerianHistorySingular',num2str(k-1),'.dat'));
    time = ( perturbedSatelliteKepler(:,1) - perturbedSatelliteKepler(1,1) )/3600;
    
    figure(k)
    for i=1:6
        
        if( i > 2 )
            multiplier = 180/pi;
        else
            multiplier = 1;
        end
        subplot(2,3,i)
        plot(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
        xlim([0 24])
        xlabel('t [hr]','interpreter','latex')
        grid on
        if( i == 1)
            ylabel('$a$ [m]','interpreter','latex')
        elseif( i == 2)
            ylabel('$e$ [-]','interpreter','latex')
        elseif( i == 3)
            ylabel('$i$ [deg]','interpreter','latex')
        elseif( i == 4)
            ylabel('$\omega$ [deg]','interpreter','latex')
        elseif( i == 5)
            ylabel('$\Omega$ [deg]','interpreter','latex')
        elseif( i == 6)
            ylabel('$\theta$ [deg]','interpreter','latex')
        end
    end
    
    if( k == 1 )
        suptitle('Earth orbit - circular')
    elseif( k == 2 )
        suptitle('Earth orbit - equatorial')
    elseif( k == 3 )
        suptitle('Earth orbit - circular and equatorial')
    end
    
    set( figure(k), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
    set( figure(k),'PaperUnits','centimeters','PaperPosition',[0 0 60 80/3]);
    saveas(gcf,strcat('singularOrbitKepler',num2str(k)),'epsc');
    
    
end

%%
for k=1:4
    perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationMeeHistorySingular',num2str(k-1),'.dat'));
    time = ( perturbedSatelliteKepler(:,1) - perturbedSatelliteKepler(1,1) )/3600;
    
    figure(k+3)
    for i=1:6
        if( i > 5 )
            multiplier = 180/pi;
        else
            multiplier = 1;
        end
        subplot(2,3,i)
        plot(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
        xlim([0 24])
        xlabel('t [hr]','interpreter','latex')
        grid on
        if( i == 1)
            ylabel('$p$ [m]','interpreter','latex')
        elseif( i == 2)
            ylabel('$f$ [-]','interpreter','latex')
        elseif( i == 3)
            ylabel('$g$ [-]','interpreter','latex')
        elseif( i == 4)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$h$ [-]','interpreter','latex')
        elseif( i == 5)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$k$ [-]','interpreter','latex')
        elseif( i == 6)
            ylabel('$L$ [deg]','interpreter','latex')
        end
        
    end
    
    if( k == 1 )
        suptitle('Earth orbit - circular')
    elseif( k == 2 )
        suptitle('Earth orbit - equatorial')
    elseif( k == 3 )
        suptitle('Earth orbit - circular and equatorial')
    elseif( k == 4 )
        suptitle('Earth orbit - retrograde equatorial')

    end
    
    
    set( figure(k+3), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
    set( figure(k+3),'PaperUnits','centimeters','PaperPosition',[0 0 60 80/3]);
    saveas(gcf,strcat('singularOrbitMee',num2str(k+3)),'epsc')
    
end

%%
for k=1:4
    perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationMeeFromKeplerHistorySingular',num2str(k-1),'.dat'));
    time = ( perturbedSatelliteKepler(:,1) - perturbedSatelliteKepler(1,1) )/3600;
    
    figure(k+7)
    for i=1:6
        if( i > 5 )
            multiplier = 180/pi;
        else
            multiplier = 1;
        end
        subplot(2,3,i)
        plot(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
        xlim([0 24])
        xlabel('t [hr]','interpreter','latex')
        grid on
        if( i == 1)
            ylabel('$p$ [m]','interpreter','latex')
        elseif( i == 2)
            ylabel('$f$ [-]','interpreter','latex')
        elseif( i == 3)
            ylabel('$g$ [-]','interpreter','latex')
        elseif( i == 4)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$h$ [-]','interpreter','latex')
        elseif( i == 5)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$k$ [-]','interpreter','latex')
        elseif( i == 6)
            ylabel('$L$ [deg]','interpreter','latex')
        end
        
    end
    
    if( k == 1 )
        suptitle('Earth orbit - circular')
    elseif( k == 2 )
        suptitle('Earth orbit - equatorial')
    elseif( k == 3 )
        suptitle('Earth orbit - circular and equatorial')
    elseif( k == 4 )
        suptitle('Earth orbit - retrograde equatorial')

    end
    
    
    set( figure(k+3), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
    set( figure(k+3),'PaperUnits','centimeters','PaperPosition',[0 0 60 80/3]);
    saveas(gcf,strcat('singularOrbitMee',num2str(k+3)),'epsc')
    
end
%%
close all
for k=1:4
    perturbedSatelliteKepler = load(strcat(dataFolder,'singlePerturbedSatellitePropagationUsmHistorySingular',num2str(k-1),'.dat'));
    time = ( perturbedSatelliteKepler(:,1) - perturbedSatelliteKepler(1,1) )/3600;
    
    figure(k+11)
    for i=1:7

        subplot(2,4,i)
        scatter(time,multiplier*perturbedSatelliteKepler(:,i+1),'LineWidth',2)
        xlim([0 24])
        xlabel('t [hr]','interpreter','latex')
        grid on
        if( i == 1)
            ylabel('$p$ [m]','interpreter','latex')
        elseif( i == 2)
            ylabel('$f$ [-]','interpreter','latex')
        elseif( i == 3)
            ylabel('$g$ [-]','interpreter','latex')
        elseif( i == 4)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$h$ [-]','interpreter','latex')
        elseif( i == 5)
            if( k == 4 )
                ylim([-1E6 1E6])
            end
            ylabel('$k$ [-]','interpreter','latex')
        elseif( i == 6)
            ylabel('$L$ [deg]','interpreter','latex')
        end
        
    end
    
    if( k == 1 )
        suptitle('Earth orbit - circular')
    elseif( k == 2 )
        suptitle('Earth orbit - equatorial')
    elseif( k == 3 )
        suptitle('Earth orbit - circular and equatorial')
    elseif( k == 4 )
        suptitle('Earth orbit - retrograde equatorial')

    end
    
    
    set( figure(k+11), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
    set( figure(k+11),'PaperUnits','centimeters','PaperPosition',[0 0 60 80/3]);
    saveas(gcf,strcat('singularOrbitMee',num2str(k+3)),'epsc')
    
end