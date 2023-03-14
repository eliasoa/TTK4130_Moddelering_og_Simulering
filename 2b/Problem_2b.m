close all
clear
clc
%% Simulation parameters
tspan = [0 10];
alpha = 10^(-3);
epsilon = [10^(-3) 10^(-6)];
y0 = [1 1 1 1];

%% System as ODE
[t_ode_e3,y_ode_e3] = ode15s(@(t,x) ODE(t,x,alpha,epsilon(1)),tspan,y0);

[t_ode_e6,y_ode_e6] = ode15s(@(t,x) ODE(t,x,alpha,epsilon(2)),tspan,y0);

%% System as DAE
M = diag([1 1 0 0]);
options = odeset('Mass',M);
[t_dae,y_dae] = ode15s(@(t,x) DAE(t,x,alpha),tspan,y0,options);


%% Plots
figure(1)
subplot(4,1,1)

hold on 
plot(t_ode_e3,y_ode_e3(:,1));
plot(t_dae,y_dae(:,1));
hold off

grid
title('$\epsilon = 10^{-3}$',Interpreter='latex')
ylabel('$x_1$',Interpreter='latex')
xlabel('t')
legend('ODE','DAE')

subplot(4,1,2)
hold on
plot(t_ode_e3,y_ode_e3(:,2));
plot(t_dae,y_dae(:,2));
hold off
grid
ylabel('$x_2$',Interpreter='latex')

subplot(4,1,3)
hold on
plot(t_ode_e3,y_ode_e3(:,3));
plot(t_dae,y_dae(:,3));
hold off
grid
ylabel('$z_1$',Interpreter='latex')

subplot(4,1,4)
hold on
plot(t_ode_e3,y_ode_e3(:,4));
plot(t_dae,y_dae(:,4));
hold off
grid
ylabel('$z_2$',Interpreter='latex')
xlabel('t')

figure(2)
subplot(4,1,1)

hold on 
plot(t_ode_e6,y_ode_e6(:,1));
plot(t_dae,y_dae(:,1));
hold off

grid
title('$\epsilon = 10^{-6}$',Interpreter='latex')
ylabel('$x_1$',Interpreter='latex')
xlabel('t')
legend('ODE','DAE')

subplot(4,1,2)
hold on
plot(t_ode_e6,y_ode_e6(:,2));
plot(t_dae,y_dae(:,2));
hold off
grid
ylabel('$x_2$',Interpreter='latex')

subplot(4,1,3)
hold on
plot(t_ode_e6,y_ode_e6(:,3));
plot(t_dae,y_dae(:,3));
hold off
grid
ylabel('$z_1$',Interpreter='latex')

subplot(4,1,4)
hold on
plot(t_ode_e6,y_ode_e6(:,4));
plot(t_dae,y_dae(:,4));
hold off
grid
ylabel('$z_2$',Interpreter='latex')
xlabel('t')
%% Functions for ode15s
function ode = ODE(t,x,alpha_ODE,epsilon_ode)

% z = x(3) & x(4)
ode = [-x(1) - x(2) - x(3);
       -x(2) - x(4);
       (1/epsilon_ode)*(0.1*x(1) - (x(1)^2+alpha_ODE)*x(3) - x(2)*x(4));
       (1/epsilon_ode)*(0.1*x(2) - (x(2)^2+alpha_ODE)*x(4))];
end
function dae = DAE(t,x,alpha_DAE)
dae = [-x(1) - x(2) - x(3);
       -x(2) - x(4);
       0.1*x(1) - (x(1)^2+alpha_DAE)*x(3) - x(2)*x(4);
       0.1*x(2) - (x(2)^2+alpha_DAE)*x(4)];
end