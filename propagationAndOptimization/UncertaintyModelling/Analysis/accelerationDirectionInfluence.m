
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontSize',12)

clc
clear all
close all

dataFolder = '../../SimulationOutput/UncertaintyModelling/';

for i=1:10
    orbitsCircular{i} = load(strcat(dataFolder,'accelerationDirectionInfluence_',num2str(i-1),'_0.dat'));
    if( i > 1 )
        differencesCircular{i-1} = orbitsCircular{i} - orbitsCircular{1};
    end
end

for i=1:10
    orbits{i} = load(strcat(dataFolder,'accelerationDirectionInfluence_',num2str(i-1),'_1.dat'));
    if( i > 1 )
        differences{i-1} = orbits{i} - orbits{1};
    end
end

%%
close all

caseLabels{1} = ' constant';
caseLabels{2} = ' $\sin(\theta)$';
caseLabels{3} = ' $\cos(\theta)$';


for testCase = 0:2    
    figure(testCase + 1)
    for i=1:3
        subplot(1,3,i)
        
        if( i == 1 )
            legend('x','y','z')
        end
        plot(orbits{1}(:,1)/86400,sqrt(differences{i+3*testCase}(:,2).^2+differences{i+3*testCase}(:,3).^2+differences{i+3*testCase}(:,4).^2))
        grid on
        
        xlim([0 14])
        
        xlabel('Time [days]')
        ylabel('$|\Delta \mathbf{r}|$ [m]')
        
        if( i == 1 )
            title( 'Radial acc.' )
        elseif( i == 2 )
            title( 'Along-track acc.' )
        elseif( i == 3 )
            title( 'Cross-track acc.' )
        end
        
    end
    
    suptitle(strcat('Impact of ',caseLabels{testCase+1},' 10$^{-8}$ m/s$^{2}$ acceleration, e=0.05'))
    
    set( figure(testCase + 1), 'Units', 'normalized', 'Position', [0,0,0.5,0.5]);
    set( figure(testCase + 1),'PaperUnits','centimeters','PaperPosition',[0 0 30 15]);
    pause(1.0)
    
    saveas(gcf,strcat('accelerationDifferenceInfluence_0_',num2str(testCase)),'png')
end


for testCase = 0:2    
    figure(testCase + 4)
    for i=1:3
        subplot(1,3,i)
        
        if( i == 1 )
            legend('x','y','z')
        end
        plot(orbitsCircular{1}(:,1)/86400,sqrt(differencesCircular{i+3*testCase}(:,2).^2+differencesCircular{i+3*testCase}(:,3).^2+differencesCircular{i+3*testCase}(:,4).^2))
        grid on
        
        xlim([0 14])
        
        xlabel('Time [days]')
        ylabel('$|\Delta \mathbf{r}|$ [m]')
        
        if( i == 1 )
            title( 'Radial acc.' )
        elseif( i == 2 )
            title( 'Along-track acc.' )
        elseif( i == 3 )
            title( 'Cross-track acc.' )
        end
        
    end
    
    suptitle(strcat('Impact of ',caseLabels{testCase+1},' 10$^{-8}$ m/s$^{2}$ acceleration, e=0.01'))
    
    set( figure(testCase + 4), 'Units', 'normalized', 'Position', [0,0,0.5,0.5]);
    set( figure(testCase + 4),'PaperUnits','centimeters','PaperPosition',[0 0 30 15]);
    pause(1.0)
    
    saveas(gcf,strcat('accelerationDifferenceInfluence_1_',num2str(testCase)),'png')
end







