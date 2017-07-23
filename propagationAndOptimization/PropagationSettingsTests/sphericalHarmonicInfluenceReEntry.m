
clc
clear all
close all

currentPath = '/home/dominic/Software/tudatBundle/tudatApplications/propagationAndOptimization/';
dataPath = currentPath;

orbit= cell(5,1);

figure(1)
for i=1:6
    orbit{i,1}=load(strcat(dataPath,'stateReEntrySphericalHarmonicCases_',num2str(i-1),'_0_0.dat'));
    %if( i == 2 )
    if( i > 1 )
        subplot(1,5,i-1)
        plot(orbit{i,1}(1:200,1),orbit{i,1}(1:200,2:4)-orbit{i-1,1}(1:200,2:4))
        
        if( i == 2 )
            title({'D/O: 16/16 w.r.t. 4/4';' '})
        elseif( i == 3 )
            title({'D/O: 4/4 w.r.t. 4/0';' '})
        elseif( i == 4 )
            title({'D/O: 4/0 w.r.t. 3/0';' '})
        elseif( i == 5 )
            title({'D/O: 3/0 w.r.t. 2/0';' '})
        elseif( i == 6 )
            title({'D/O: 2/0 w.r.t. 0/0';' '})
        end
        
        xlabel('Time [s]')
        if( i== 2 )
            ylabel('Position [m]')
            legend('x','y','z','Location','NorthWest')
        end
        grid on
    end
    
    
    
    
end
suptitle({'Entry capsule position difference for different maximum D/O of Earth gravity field';''})

set(figure(1), 'Units', 'normalized', 'Position', [0,0,0.66,0.5]);
set(figure(1),'PaperUnits','centimeters','PaperPosition',[0 0 40 20 ])
set(figure(1),'PaperPositionMode','auto')
