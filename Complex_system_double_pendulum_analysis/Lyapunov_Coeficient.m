%% Introduction: double pendulum
%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;
%% Coeficient Lyapunov

tspan = linspace(0, 10, 1000);
angle = linspace(1, 60, 30);
result = zeros(length(tspan), 2, length(angle)); % Store distances

for i = 1:length(angle)
    theta0 = angle(i) * pi / 180;
    initial_condition = [theta0; 0; 0; 0];
    
    [~, z] = ode45(@(t, theta) equations(t, theta, L1, L2, g, m1, m2), tspan, initial_condition);
    [~, z_aprox] = ode45(@(t, theta) equations_aprox(t, theta, L1, L2, g, m1, m2), tspan, initial_condition);
    
    % Interpolate if needed
    z_aprox_interp = interp1(tspan, z_aprox, tspan);
    
    result(:,:,i) = compare(z(:,1), z(:,3), z_aprox_interp(:,1), z_aprox_interp(:,3));
end
%% PLots

% Compute mean error across time and both pendulums
mean_error = squeeze(mean(mean(result, 2), 1)); % [1 x numAngles]

% Plot
figure(1);
plot(angle, mean_error, 'r-o', 'LineWidth', 2);
xlabel('Initial Angle (degrees)');
ylabel('Mean Absolute Error');
title('Accuracy of Linearized vs Nonlinear Double Pendulum');
grid on;

figure(2);

% Assuming result is a 3D matrix: [timeSteps x 2 x numAngles]
% Let's compute two metrics for comparison:
mean_error1 = squeeze(mean(result(:,1,:), 1)); % Error for pendulum 1
mean_error2 = squeeze(mean(result(:,2,:), 1)); % Error for pendulum 2

% Plot both error curves
plot(angle, mean_error1, 'b-o', 'LineWidth', 2); 
hold on;
plot(angle, mean_error2, 'g-o', 'LineWidth', 2);

xlabel('Initial Angle (degrees)');
ylabel('Mean Absolute Error');
title('Comparison of Linearized vs Nonlinear Models');
legend('Pendulum 1', 'Pendulum 2');
grid on;

%% Comparetion
function distance = compare(z1, z2, z_aprox1, z_aprox2)
    distance = [abs(z1 - z_aprox1), abs(z2 - z_aprox2)];

end
%% Equation
function dy = equations(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    dtheta1 = theta(2);
    theta2 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1^2, m2 * L1 * L2 * cos(theta1 - theta2);
        m2 * L1 * L2 * cos(theta1 - theta2), m2 * L2^2];

    B = [-m2 * L1 * L2 * sin(theta1 - theta2) * dtheta2^2 - (m1 + m2) * g * L1 * sin(theta1);
     m2 * L1 * L2 * sin(theta1 - theta2) * dtheta1^2 + m2 * g * L2 * sin(theta2)];

    sol = A \ B;

    dtheta1_sol = sol(1);
    dtheta2_sol = sol(2);
    
    
    dy = zeros(4,1);
    dy(1) = dtheta1;
    dy(2) = dtheta1_sol;
    dy(3) = dtheta2;
    dy(4) = dtheta2_sol;

end
%% Equation_aprox
function dy_aprox = equations_aprox(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    dtheta1 = theta(2);
    theta2 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1, m2 * L2;
         m2 * L1,        m2 * L2];

    B = [- (m1 + m2) * g * theta1;
         - m2 * g * theta2];

    sol = A \ B;

    dtheta1_aprox = sol(1);
    dtheta2_aprox = sol(2);

    dy_aprox = zeros(4,1);
    dy_aprox(1) = dtheta1;
    dy_aprox(2) = dtheta1_aprox;
    dy_aprox(3) = dtheta2;
    dy_aprox(4) = dtheta2_aprox;
end