clear
clc

alfa = pi;

C = @(alfa) [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)];

position = [2; 1];

    
figure(1)
hold on
grid on
plot(0, 0, 'xr')           % origin
plot(position(1), position(2), 'xb');

steps = 64;

for i = 1:steps - 1
    new_position = C(alfa*i/(steps/2))*position;
    plot(new_position(1), new_position(2), 'ok')
end

scale = 4;
axis equal
axis([-scale scale -scale scale])
legend('Origin','Starting point')

