
%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;

%Condition
t_init = [0 10];
angle = linspace(1,45,44);

%% Initial angle variation for the double pendulum

E_total_average = zeros(length(angle),2);

for i=1:length(angle)
    initial_condition_aprox = [angle(i) * pi / 180; 0; 0; 0];
    [t, z] = ode45(@(t, theta)equations(t, theta, L1, L2, g, m1, m2), t_init, initial_condition_aprox, odeset);
    
    Ek = zeros(length(t), 1);
    Ep = zeros(length(t), 1);
    E_total= zeros(length(t),1);
    for j = 1:length(t)
    Ek(j) = 0.5*(m1 + m2)*L1^2*z(j,2)^2 + 0.5*m2*L2^2*z(j,4)^2 + m2*L1*L2*z(j,2)*z(j,4)*cos(z(j,1)-z(j,3));
    Ep(j) = -(m1 + m2)*g*L1*cos(z(j,1)) - m2*g*L2*cos(z(j,3));
    E_total(j) = Ek(j) + Ep(j);
    end
    
    E_total_average(i,1) = max(E_total);
    E_total_average(i,2) = min(E_total);
end 

figure(1)
plot(angle, E_total_average(:,1), 'b-')
hold on
plot(angle,E_total_average(:,2), 'r-')
xlabel('Initial angle (degrees)'); ylabel('Average Total Energy');
legend('Max Energy', 'Min Energy')
grid on

energy_range = E_total_average(:,1) - E_total_average(:,2);
figure(2)
plot(angle, energy_range, 'k-', 'LineWidth', 2)
xlabel('Initial angle (degrees)');
ylabel('Energy Range (J)');
title('Energy Spread vs Initial Angle');
grid on

%% Initial angle variation for the double pendulum (approximation)

E_total_average_approx = zeros(length(angle),2);

for i=1:length(angle)
    initial_condition_aprox = [angle(i) * pi / 180; 0; 0; 0];
    [t_aprox, z_aprox] = ode45(@(t, theta)equations_aprox(t, theta, L1, L2, g, m1, m2), t_init, initial_condition_aprox, odeset);
    
    Ek_aprox = zeros(length(t_aprox), 1);
    Ep_aprox = zeros(length(t_aprox), 1);
    E_total_aprox = zeros(length(t_aprox),1);
    for j = 1:length(t_aprox)
    Ek_aprox(j) = 0.5*(m1 + m2)*L1^2*z_aprox(j,2)^2 + 0.5*m2*L2^2*z_aprox(j,4)^2 + m2*L1*L2*z_aprox(j,2)*z_aprox(j,4)*cos(z_aprox(j,1)-z_aprox(j,3));
    Ep_aprox(j) = -(m1 + m2)*g*L1*cos(z_aprox(j,1)) - m2*g*L2*cos(z_aprox(j,3));
    E_total_aprox(j) = Ek_aprox(j) + Ep_aprox(j);
    end
    
    E_total_average_approx(i,1) = max(E_total_aprox);
    E_total_average_approx(i,2) = min(E_total_aprox);
end 

figure(3)
plot(angle, E_total_average_approx(:,1), 'b-')
hold on
plot(angle,E_total_average_approx(:,2), 'r-')
xlabel('Initial angle (degrees)'); ylabel('Average Total Energy (approx)');
legend('Max Energy', 'Min Energy')
grid on

energy_range = E_total_average_approx(:,1) - E_total_average_approx(:,2);
figure(4)
plot(angle, energy_range, 'k-', 'LineWidth', 2)
xlabel('Initial angle (degrees)');
ylabel('Energy Range (J)');
title('Energy Spread vs Initial Angle (approx)');
grid on
%% Function
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