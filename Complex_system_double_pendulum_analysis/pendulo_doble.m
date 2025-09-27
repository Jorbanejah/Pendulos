%%Introduction: double pendulum
%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;

%Condition
t_init = [0 10];
initial_condition = [45 * pi / 90; 0; 0; 0];
initial_condition_aprox = [45 * pi / 90; 0; 0; 0];
[t, z] = ode45(@(t, theta)equations(t, theta, L1, L2, g, m1, m2), t_init, initial_condition);
[t_aprox, z_aprox] = ode45(@(t, theta)equations_aprox(t, theta, L1, L2, g, m1, m2), t_init, initial_condition_aprox);

%%
%Animation: with a little tricky


if (length(t) >= length(t_aprox))
    nFrames = length(t_aprox);
else
    nFrames = length(t);
end

for i = 1:nFrames
    clf;
    x = L1 * sin(z(i,1));
    y = - L1 * cos(z(i,1));
    x1 = x + L2 * sin(z(i,3));
    y1 = y - L2 * cos(z(i,3));

    x_aprox = L1 * sin(z_aprox(i,1));
    y_aprox = - L1 * cos(z_aprox(i,1));
    x1_aprox = x_aprox + L2 * sin(z_aprox(i,3));
    y1_aprox = y_aprox - L2 * cos(z_aprox(i,3));

    hold on

    % Plot exact pendulum
    plot([0 x], [0 y], 'k-', 'LineWidth', 2);
    plot([x x1], [y y1], 'k-', 'LineWidth', 2);
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
    plot(x1, y1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','g');

    % Plot approximated pendulum (shifted right)
    plot([2 2+x_aprox], [0 y_aprox], 'k-', 'LineWidth', 2);
    plot([x_aprox+2 x1_aprox+2], [y_aprox y1_aprox], 'k-', 'LineWidth', 2);
    plot(x_aprox+2, y_aprox, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
    plot(x1_aprox+2, y1_aprox, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','g');

    % Ground line
    plot([-3 5], [0 0], 'k', 'LineWidth', 2);

    axis equal; axis([-2, 4, -3, 1]);
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
%%
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