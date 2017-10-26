clc
clear all
close all

folder = '../../SimulationOutput/EnvironmentModels/';


orbit= cell(3,3);
orbit_full= cell(3,3);

figure(1)
for i=1:3
    for j=1:4
        orbit{i,j}=load(strcat(folder,'stateMoonOrbiterRotationCases_full_',num2str(i-1),'_',num2str(j-1),'.dat'));
        
        if( j > 1 )
            subplot(3,3,(i-1)*3+(j-1))
            plot(orbit{i,j}(:,1),orbit{i,j}(:,2:4)-orbit{i,1}(:,2:4))
        end
        
    end
end
%%
close all
figure(2)
for i=1:3
    for j=1:4
        %orbit{i,j}=load(strcat(folder,'stateMoonOrbiterRotationCases_full_',num2str(i-1),'_',num2str(j-1),'.dat'));
        
        if( i > 1 )
            subplot(2,4,(i-2)*4+(j))
            plot(orbit{i,j}(:,1)/86400,orbit{i,j}(:,2:4)-orbit{1,j}(:,2:4))
            if( i == 2 )
                ylim([-150 150])
            elseif( i == 2 )
                ylim([-5 5])
            end
            grid on
            if( i == 3 )
                xlabel('Time [days]')
            end
            if( j == 1 )
                ylabel('Position difference [m]')
            end
                                  
            if( i== 2 )
                if( j == 1 )
                    title('D/O=16/16')
                elseif( j == 2 )
                    title('D/O=4/4')                    
                elseif( j == 3 )
                    title('D/O=2/2')                    
                elseif( j == 4 )
                    title('D/O=2/0')                    
                end
            end
        end
    end
end

for figureIndex = 2
    set(figure(figureIndex), 'Units', 'normalized', 'Position', [0,0,0.75,0.66]);
    set(figure(figureIndex),'PaperUnits','centimeters','PaperPosition',[0 0 45 22.5 ])
    set(figure(figureIndex),'PaperPositionMode','auto')
end