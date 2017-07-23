

clc
clear all
close all

currentPath = '/home/dominic/Software/tudatBundle/tudatApplications/propagationAndOptimization/';
dataPath = strcat(currentPath,'EphemerisSettings/');

orbit= cell(2,6);

figure(1)
for i=1:12
    column = i;
    row = 1;
    if( column > 6 )
        column = column - 6;
        row = row + 1;
    end
    orbit{row,column}=load(strcat(dataPath,'stateMoonOrbiterSphericalHarmonicCases_full_',num2str(i-1),'.dat'));
    
    
    if( column > 1 )
        subplot(2,5,(row - 1)*5+column-1)
        plot(orbit{row,column}(:,1),orbit{row,column}(:,2:4)-orbit{row,column-1}(:,2:4))
        if( row == 1 )
            if( column == 2 )
                title({'D/O: 128/128 w.r.t. 32/32';' '})
            elseif( column == 3 )
                title({'D/O: 32/32 w.r.t. 8/8';' '})
            elseif( column == 4 )
                title({'D/O: 8/8 w.r.t. 2/2';' '})
            elseif( column == 5 )
                title({'D/O: 2/2 w.r.t. 2/0';' '})
            elseif( column == 6 )
                title({'D/O: 2/0 w.r.t. 0/0';' '})
            end
        end
    end
    
    
    
    if( row == 2 )
        xlabel('Time [s]')
    end
    if( column == 2 )
        ylabel('Position difference [m]')
    end
end
suptitle('Lunar orbiter position difference for different maximum D/O of lunar gravity field')
%%
set(figure(1), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
set(figure(1),'PaperUnits','centimeters','PaperPosition',[0 0 60 22.5 ])
set(figure(1),'PaperPositionMode','auto')
%%
figure(2)
for i=1:12
    column = i;
    row = 1;
    if( column > 6 )
        column = column - 6;
        row = row + 1;
    end
    
    if( column < 6 )
        subplot(2,6,(row - 1)*6+column)
        plot(orbit{row,column}(:,1),orbit{row,column}(:,2:4)-orbit{row,6}(:,2:4))
    end
end