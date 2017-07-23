clc
close all
clear all

load('/home/dominic/Downloads/outputMap.dat')

load('/home/dominic/Downloads/inertialTorqueMap.dat')
load('/home/dominic/Downloads/inertialAngularMomentumMap.dat')
load('/home/dominic/Downloads/inertialAngularMomentumDerivativeMap.dat')
load('/home/dominic/Downloads/inertialAngularMomentumDerivativeErrorMap.dat')

nRows = size(inertialTorqueMap, 1);
times = inertialAngularMomentumMap(:,1);

figure(4)
subplot(1,3,1)
plot(inertialTorqueMap(:,1),inertialTorqueMap(:,2:4),'-x')

for i=1:3
    angularMomentumRate(:,i) = (inertialAngularMomentumMap(2:nRows,i+1) - inertialAngularMomentumMap(1:(nRows-1),i+1))./(times(2:nRows) - times(1:(nRows-1)));
end
grid on

subplot(1,3,2)
plot(times(2:nRows),angularMomentumRate,'-x')
grid on

subplot(1,3,3)
plot(times(2:nRows),(times(2:nRows) - times(1:(nRows-1))),'-x')
grid on

figure(5)
subplot(1,2,1)
plot(inertialAngularMomentumDerivativeErrorMap(:,1),inertialAngularMomentumDerivativeErrorMap(:,2:4),'-x')
grid on

subplot(1,2,2)
plot(inertialAngularMomentumDerivativeErrorMap(:,1),sqrt(inertialAngularMomentumDerivativeErrorMap(:,2).^2+inertialAngularMomentumDerivativeErrorMap(:,3).^2+inertialAngularMomentumDerivativeErrorMap(:,4).^2),'-x')
grid on

%%
times =  outputMap(:,1);
torques = outputMap(:,8:10);
inertiaTensor = zeros(3,3);
inertiaTensor(1,1) = outputMap(1,11);
inertiaTensor(2,2) = outputMap(1,12);
inertiaTensor(3,3) = outputMap(1,13);
angularVelocity = outputMap(:,14:16);

nRows = size(outputMap, 1);
for i=1:3
    angularVelocityRate(:,i) = inertiaTensor(i,i) * (angularVelocity(2:nRows,i) - angularVelocity(1:(nRows-1),i))./(times(2:nRows) - times(1:(nRows-1)));
end

rhs = zeros(nRows,3);

for i=1:nRows
    rhs1(i,1:3)= ((-cross(angularVelocity(i,:)',inertiaTensor*(angularVelocity(i,:)'))))';
    rhs2(i,1:3)= ((torques(i,:)'))';
    rhsTot= rhs1 + rhs2;
end

figure(1)

subplot(1,4,1)
plot(outputMap(:,1),rhs1)
subplot(1,4,2)
plot(outputMap(:,1),rhs2)
subplot(1,4,3)
plot(outputMap(:,1),rhsTot)
subplot(1,4,4)
plot(times(2:nRows),angularVelocityRate)
