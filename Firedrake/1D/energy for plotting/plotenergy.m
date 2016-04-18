clc
close all


T = dlmread('energy_100_2_16_16.txt');
figure 
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')