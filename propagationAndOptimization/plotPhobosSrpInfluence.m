clc
close all
clear all

%%

figure(1)
dataFolder = 'phobosOPropagationResults/';
perturbedPhobos = load(strcat(dataFolder,'phobosPropagationHistorySingular0.dat'));
unperturbedPhobos = load(strcat(dataFolder,'phobosPropagationHistorySingular1.dat'));
orbitDifference = perturbedPhobos - unperturbedPhobos;

time = ( perturbedPhobos(:,1) - unperturbedPhobos(1,1) )/86400;

plot(time,orbitDifference(:,2),'LineWidth',2)
hold on
plot(time,orbitDifference(:,3),'LineWidth',2)
plot(time,orbitDifference(:,4),'LineWidth',2)
grid on

xlabel('$t$ [s]','interpreter','latex')
ylabel('$\Delta r$ [m]','interpreter','latex')
legend('x','y','z')

xlim([0 10*365.25])
set( gcf, 'Units', 'normalized', 'Position', [0,0,0.5,0.75]);
set( gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 30]);
pause(2.0)
saveas(gcf,strcat('phobosSrpInfluence'),'epsc');
%%
perturbedKeplerPhobos = load(strcat(dataFolder,'phobosPropagationKeplerianHistorySingular0.dat'));
unperturbedKeplerPhobos = load(strcat(dataFolder,'phobosPropagationKeplerianHistorySingular1.dat'));
orbitDifferenceKepler = perturbedKeplerPhobos - unperturbedKeplerPhobos;

figure(2)

for i=1:6
    if( i > 2 )
        multiplier = 180/pi;
    else
        multiplier = 1;
    end
    subplot(2,3,i)
    plot(time,multiplier*orbitDifferenceKepler(:,i+1),'LineWidth',2)
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
