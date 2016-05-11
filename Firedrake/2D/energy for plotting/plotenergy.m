clc
close all


T = dlmread('fv-16-16-1000.txt');
figure 
semilogy( T(:,1), T(:,2))
set(gca,'FontSize', 14)
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')