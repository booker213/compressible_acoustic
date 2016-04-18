clc
close all


T = dlmread('fv-16-16-10.txt');
figure 
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')