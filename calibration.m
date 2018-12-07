%   Function used to calibrate the models, requires following inputs:
%   List of actual positions
%   p = [[x1, y1, z1; x2, y2, z2; ... etc]]

%   Filenames in the form of 'J_A_static_1'

function th_estimate = calibration(p)
    y = [];
    G = [];
    for i = 1:size(p,2)
        filename = sprintf('J_A_static_%i.txt',i);
        data = load_data(filename);
        y_current = data(:);

        p_current = p(:,i);
        H = (3*p_current*p_current'-norm(p_current).^2*eye(3))/norm(p_current).^5; 
        G_current = [eye(3) H];

        G = [G; repmat(G_current, size(y_current,1) / 3, 1)];
        y = [y; y_current];
    end
    % Now the estimate can be calculated:
    th_estimate = (G'*G)\G'*y;
    
    
    
    %{

y = [];
G = [];
for i = 1:size(p,2)
    filename = sprintf('sensorLog_%i.txt',i);
    data = load_data(filename);
    y_current = data(:);
    
    p_current = p(:,i);
    H = (3*p_current*p_current'-norm(p_current).^2*eye(3))/norm(p_current).^5; 
    G_current = [eye(3) H];

    G = [G; repmat(G_current, size(y_current,1) / 3, 1)];
    y = [y; y_current];
end
%}
end