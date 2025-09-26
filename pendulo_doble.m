%%Introduction: double pendulum
%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;
%angle = linspace(0, 45, 45);
%Condition
t_init = [0 2];
initial_condition = [12 * pi / 90; 0; 0; 0];
%initial_condition_aprox = [12 * pi / 90; 0; 0; 0];
[t, z] = ode45(@(t, theta)equations(t, theta, L1, L2, g, m1, m2), t_init, initial_condition);
%[t_aprox, z_aprox] = ode45(@(t, theta)equations_aprox, t_init, initial_condition_aprox);

%%
%Animation
nFrames = length(t);
%interpol = interp1(t_aprox, z_aprox, nFrames);
if (nFrames > 3000) 
    nFrames = nFrames - 1000; 
end
for i = 1:nFrames
    clf;
    x = L1 * sin(z(i,1));
    y = - L1 * cos(z(i,1));
    x1 = x + L2 * sin(z(i,3));
    y1 = y - L2 * cos(z(i,3));
    hold on
    plot([0 x], [0 y], 'k-', 'LineWidth', 2);
    plot([x x1], [y y1], 'k -', 'LineWidth', 2);
    plot([-3 3], [0 0], 'k', 'LineWidth', 2);
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'g');
    plot(x1, y1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','g', 'MarkerEdgeColor', 'g');

    axis equal; axis([-2, 2, -3, 1]);
    drawnow;
end
%%
function dy = equations(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    dtheta1 = theta(2);
    theta2 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1^2, m2 * L1 * L2 * cos(theta1 - theta2);
        m2 * L1 * L2 * cos(theta1 - theta2), m2 * L2^2];

    B = [-m2 * L1 * L2 * dtheta1^2 * sin(theta1 - theta2) - (m1 + m2) * g * L1 * sin(theta1);
        m2 * L1 * L2 * dtheta2^2 * sin(theta1 - theta2) + m2 * g * L2 * sin(theta2)];

    sol = A \ B;

    dtheta1_sol = sol(1);
    dtheta2_sol = sol(2);
    
    
    dy = zeros(4,1);
    dy(1) = dtheta1;
    dy(2) = dtheta1_sol;
    dy(3) = dtheta2;
    dy(4) = dtheta2_sol;

end
%%
function dy_aprox = equations_aprox(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    theta2 = theta(2);
    dtheta1 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1,  L2* m2;
         L1, L2];
    B = [- (m1 + m2) * g * theta1;
         - g * theta_2];

    sol = A \ B;

    dtheta1_aprox = sol(1);
    dtheta2_aprox = sol(2);

    dy_aprox = zeros(4,1);
    dy_aprox(1) = dtheta1;
    dy_aprox(2) = dtheta1_aprox;
    dy_aprox(3) = dtheta2;
    dy_aprox(4) = dtheta2_aprox;
    
end