clear

%System parameters
K = 22200;
B_E = 10;
y_s = 1;
v_E = 0.05;
d_E = 0.03;
d_F = 0.04;
d_M = 0.1;
d_s = 0.12;
v = 0.49;
eta = 0.7;

%Spacing, day amount, f(x)
h = 10^-1;
D = 50;
f = @(x, u) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            u - d_s * x(4)];

u = zeros(1, D/h + 1) + 10^5;
%Solution
x = RK4U(f, u, h, D, [21910; 5587; 13419; 0]);

T = 0:h:D;
%Plot of E, M and F
figure;
plot(T, x(1,:), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T, x(2,:), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T, x(3,:), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of E, M and F')
pbaspect([2 1 1]); 
hold off;

%Plot of M_s
figure;
plot(T, x(4,:), 'Color', [0.4660 0.6740 0.1880], 'DisplayName', ...
    'Sterilized male mosquitoes M_s', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of M_s');
pbaspect([2 1 1]);

%Plot of the control u
figure;
plot(T, u(1:D/h + 1), 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Control function u', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Release rate');
legend;
title('Plot of the control u');
pbaspect([2 1 1]); 
%%
clear

h = 10^-1;
D = 50;

N = D/h + 1;
fun = @Constraints;
x0 = zeros(4*N, 1);

X = fsolve(fun, x0);
%X = [[21910; 5587; 13419; 0], [X(1:N)'; X(N+1:2*N)'; X(2*N+1:3*N)'; X(3*N+1:4*N)']];
X = [X(1:N)'; X(N+1:2*N)'; X(2*N+1:3*N)'; X(3*N+1:4*N)'];

T=0:h:D;
%Plot of E, M and F
figure;
plot(T, X(1,:), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T, X(2,:), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T, X(3,:), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of E, M and F')
pbaspect([2 1 1]); 
hold off;

%Plot of M_s
figure;
plot(T, X(4,:), 'Color', [0.4660 0.6740 0.1880], 'DisplayName', ...
    'Sterilized male mosquitoes M_s', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of M_s');
pbaspect([2 1 1]);