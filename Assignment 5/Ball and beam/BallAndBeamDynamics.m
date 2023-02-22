function x_dot = BallAndBeamDynamics(t, state, parameters)

    q = state(1:2);
    dq = state(3:4);

    T = 200*(q(1)-q(2)) + 70*(dq(1)-dq(2));

    [W, RHS] = BallAndBeamODEMatrices(state,T,parameters);
    
    x_dot = [dq;W\RHS];
end