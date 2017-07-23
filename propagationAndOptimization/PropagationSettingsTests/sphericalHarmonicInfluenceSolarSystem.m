clc
clear all
close all

currentPath = '/home/dominic/Software/tudatBundle/tudatApplications/propagationAndOptimization/';
dataPath = strcat(currentPath,'EphemerisSettings/');

orbit= cell(2,1);

figure(1)
for i=1:2
    orbit{i,1}=load(strcat(dataPath,'stateSolarSystemSphericalHarmonicCases_full_',num2str(i-1),'.dat'));
    if( i == 2 )
        for j=1:7
            subplot(1,7,j)
            plot(orbit{i,1}(:,1)/(86400*365.25),orbit{2,1}(:,(2+(j-1)*6):(2+(j-1)*6)+3)-orbit{1,1}(:,(2+(j-1)*6):(2+(j-1)*6)+3))
            if( j == 1 )
                title('Sun')
            elseif( j == 2 )
                title('Mercury')
            elseif( j == 3 )
                title('Venus')
            elseif( j == 4 )
                title('Earth')
            elseif( j == 5 )
                title('Moon (w.r.t Earth)')
            elseif( j == 6 )
                title('Mars')
            elseif( j == 7 )
                title('Jupiter')
            end
            
            
            xlabel('Time [years]')
            if( j == 1 )
                ylabel('Position difference [m]')
            end
        end
        
        
    end
end
suptitle('Influence of Solar J2 on solar system bodies')

set(figure(1), 'Units', 'normalized', 'Position', [0,0,1,0.66]);
set(figure(1),'PaperUnits','centimeters','PaperPosition',[0 0 60 22.5 ])
set(figure(1),'PaperPositionMode','auto')
