
clc
clear all
close all

currentPath = '/home/dominic/Software/tudatBundle/tudatApplications/propagationAndOptimization/';
dataPath = currentPath;

orbit= cell(5,1);

headingAngleValues = cell(3,1);
headingAngleValues{1,1}=', initial heading = 0 deg';
headingAngleValues{2,1}=', initial heading = 45 deg';
headingAngleValues{3,1}=', initial heading = 90 deg';

%for caseNumber = 2
for initialHeading = 1:3
    figure(initialHeading)
    for initialLatitude = 0:1:9
        for initialLongitude = 0:1:19
            for i=1:8
                %if( caseNumber == 1 )
                %    orbit{i,1}=load(strcat(dataPath,'stateReEntrySphericalHarmonicCases_',num2str(initialLatitude),'_',num2str(initialLongitude),'_',num2str(i-1),'.dat'));
                %elseif( caseNumber == 2 )
                orbit{i,1}=load(strcat(dataPath,'bodyFixedStateReEntrySphericalHarmonicCases_',num2str(initialLatitude),'_',num2str(initialLongitude),'_',num2str(initialHeading-1),'_',num2str(i-1),'.dat'));
                %end    %if( i == 2 )
                
                if( i > 1 )
                    
                    s1=size(orbit{i-1,1});
                    s2=size(orbit{i,1});
                    nRows = min(s1(1),s2(1));
                    subplot(1,7,i-1)
                    difference = sqrt((orbit{i,1}(1:nRows,2)-orbit{i-1,1}(1:nRows,2)).^2+(orbit{i,1}(1:nRows,3)-orbit{i-1,1}(1:nRows,3)).^2+(orbit{i,1}(1:nRows,4)-orbit{i-1,1}(1:nRows,4)).^2);
                    plot(orbit{i,1}(1:nRows,1),difference)
                    if( initialLatitude == 0 && initialLongitude == 0 )
                        if( i == 2 )
                            title({'D/O: 32/32 w.r.t. 16/16';' '})
                        elseif( i == 3 )
                            title({'D/O: 16/16 w.r.t. 8/8';' '})
                        elseif( i == 4 )
                            title({'D/O: 16/16 w.r.t. 8/8';' '})
                        elseif( i == 5 )
                            title({'D/O: 4/4 w.r.t. 4/0';' '})
                        elseif( i == 6 )
                            title({'D/O: 4/0 w.r.t. 3/0';' '})
                        elseif( i == 7 )
                            title({'D/O: 3/0 w.r.t. 2/0';' '})
                        elseif( i == 8 )
                            title({'D/O: 2/0 w.r.t. 0/0';' '})
                        end
                        
                        xlabel('Time [s]')
                        if( i== 2 )
                            ylabel('3-Dimensional Position Difference [m]')
                            %legend('x','y','z','Location','NorthWest')
                        end
                        grid on
                        hold on
                    end
                    xlim([0 300])
                end
            end
        end
    end
    suptitle({strcat('Entry capsule position difference for different maximum D/O of Earth gravity field',headingAngleValues{initialHeading});''})
    
    set(figure(initialHeading), 'Units', 'normalized', 'Position', [0,0,0.66,0.5]);
    set(figure(initialHeading),'PaperUnits','centimeters','PaperPosition',[0 0 40 20 ])
    set(figure(initialHeading),'PaperPositionMode','auto')
end
