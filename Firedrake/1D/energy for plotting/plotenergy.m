clc
close all


T = dlmread('energy_1000_periods_fv_16_16.txt');
figure 
set(gca,'fontsize', 18)
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')