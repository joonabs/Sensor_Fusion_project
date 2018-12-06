clear
clc

p_calib = [30 30 10; 0 30 10; -30 30 10; -30 0 10; -30 -30 10; 0 -30 10;
            30 -30 10; 30 0 10; 45 45 10; 60 0 10];
p_verif = [60 0 10; 30 30 10; 45 45 10; 30 45 10; 0 60 10; 45 30 10; -30 15 10; -30 -15 10];

plot_points(p_calib)
plot_points(p_verif)

function plot_points(p)
    labels = cellstr(num2str([1:size(p,1)]'));

    figure();
    hold on
    grid on
    axis([-40 70, -40, 70]);

    rectangle('Position',[-2.5 -5 5 10])

    plot(0, 0, 'xr')

    plot(p(:,1), p(:,2), 'xk');
    text(p(:,1), p(:,2), labels, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    for i = 1:size(p, 1)
        text(p(i,1),p(i,2),['(' num2str(p(i,1)) ',' num2str(p(i,2)) ')'], ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    end

    t = text(5,0,'Phone');
    s = t.FontSize;
    t.FontSize = 12;
end