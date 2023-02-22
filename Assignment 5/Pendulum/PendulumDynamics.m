function  x_dot = PendulumDynamics(t, state, parameters)   
    % [x; theta_1; theta_2]
    q = state(1:3);
    % [x_dot; theta_1_dot; theta_2_dot]
    qd = state(4:6);
    % F = -10*x - x_dot
    F = -10*q(1)-qd(1);
    
    % x_dot = q_dot; d^2(L^-1)/dq_dot^2) ( Q + dL/dq ..... )
    %                detter er W          detter er RHS       
    [W,RHS] = PendulumODEMatrices(state,F,parameters);

    x_dot = [qd;W\RHS];
end
