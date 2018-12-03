y = [];
G = [];
for i = 1:size(p,2)
    filename = sprintf('sLog_pos%i.txt',i);
    data = load_data(filename);
    y_current = data(:);
    
    p_current = p(:,i);
    H = (3*p_current*p_current'-norm(p_current).^2*eye(3))/norm(p_current).^5; 
    G_current = [eye(3) H];

    G = [G; repmat(G_current, size(y_current,1) / 3, 1)];
    y = [y; y_current];
end