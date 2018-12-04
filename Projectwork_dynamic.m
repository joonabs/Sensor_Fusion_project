%   Projectwork dynamic
%
%
%   Aleksi Mäkinen and Joonas Isometsä
clear
clc

%   Initial parameters
%   

%   The corresponding known actual positions:
%   [x1, y1, z1; x2, y2, z2; ... etc]
p = []';

sigma_v = 0.01;
sigma_alpha = 0,1;

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

% Now the estimate can be calculated:
th_estimate = (G'*G)\G'*y;
% The estimate is now precalculated m_T !
m_T = th_estimate;



%   Solved symbolic jacobian:
function G = jacobian_symbolic()
    syms p_x p_y p_z m_T_1 m_T_2 m_T_3 alpha real;
    
    p = [p_x ; p_y ; p_z];
    m_T = [m_T_1; m_T_2; m_T_3];
    C = [cos(alpha) -sin(alpha) 0;...
        sin(alpha) cos(alpha) 0;...
        0 0 1];
    
    y = ((3*p*p' - norm(p).^2*eye(3))/norm(p).^5)*C*m_T
    G = jacobian(y, [p_x, p_y]);
    
end
%   Function to insert and calculate numeric jacobian
function G = jacobian_substitute(G, x, m_T)
    p_x = x(1);
    p_y = x(2);
    p_z = 10;
    m_T_1 = m_T(4);
    m_T_2 = m_T(5);
    m_T_3 = m_T(6);
    alpha = 
    
    G = subs(G);
    G = double(G);
end
