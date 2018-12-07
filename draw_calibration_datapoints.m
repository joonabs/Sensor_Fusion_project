%   Used to draw measurements
clc
clear all
close all
p = [30, 0, 10; 
    15, 25, 10;
    10, 45, 10;
    -30, 45, 10;
    -30, 0, 10;
    -20, -30, 10;
    15, -25, 10]';

figure()
hold on
grid on
axis([-40,40, -50, 50])
plot(0,0,'xk')
text(0,0,'(0,0)')
plot(p(1,:), p(2,:), 'rx')
names = {'1 (30,0,10)';'2 (15,25,10)';'3 (10,45,10)';'4 (-30,45,10)';'5 (-30,0,10)';'6 (-20,30,10)';'7 (15, -25, 10)'}; 
%names ={'x1';'x2';'x3'}
title('The measurement points used in calibration')
xlabel('[cm]');
ylabel('[cm]');
text(p(1,:),p(2,:),names,'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

