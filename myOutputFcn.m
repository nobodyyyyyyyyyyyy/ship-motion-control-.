function myOutputFcn(t, y, ~)
    global Flag;  % Declare the global variable
    L = 7;
    
    % Extract x and y
    u = y(1);
    v = y(2);
    r = y(3);

    % Calculate u*v*r at each time step
    V = sqrt(v^2 + u^2);
    beta = atan2(v, u);
    beta_p = beta - (-0.5 * r * L / V);

    % Check if beta_p is greater than or equal to 0 and update the global variable
    if beta_p >= 0
        Flag = 2;  % Set the global variable to 2
    end
end
