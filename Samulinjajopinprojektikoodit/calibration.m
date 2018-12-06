function [th_estimate, R_estimate] = calibration()
    %calibration_points = [20 0 11; 40 0 11; 50 0 11; 0 30 11; 30 30 11]';
    calibration_points = [-20 0 11; -30 0 11; -40 0 11; -20 -30 11; 20 -30 11]';
    
    
    %calibration_points = [-20 0 11; -30 0 11; -40 0 11; -20 -30 11]';
    
    %calibration_points_J_A = [30 0 11; 15 25 11; 0 45 11; -30 45 11; -30 0 11; -20 -30 11; 15 -25 11; 0 -55 11];
    %calibration_points = calibration_points_J_A;
    
    y = [];
    data_all = [];
    G = [];
    %for i = 2:size(calibration_points,2)
    for i = 1:size(calibration_points,2)
        %filename = sprintf('Calibration_data_%i.txt',i);
        
        filename = sprintf('Verification_data_%i.txt',i);
        %filename = sprintf(,i);
        data = load_data(filename);
        %data(n:n:end) = [];
        y_current = data(:);

        p_current = calibration_points(:,i);
        H = (3*p_current*p_current'-norm(p_current).^2*eye(3))/norm(p_current).^5; 
        G_current = [eye(3) H];

        G = [G; repmat(G_current, size(y_current,1) / 3, 1)];
        y = [y; y_current];
        data_all = [data_all data];
    end
    
    R_estimate = kron(eye(size(data_all,2)), cov(data_all'));
    th_estimate = (G'*inv(R_estimate)*G)\(G'*inv(R_estimate)*y)
end