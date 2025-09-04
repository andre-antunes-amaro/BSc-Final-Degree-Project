%Modelling of wild mosquito population dynamics with allee effect%

%Initial condition: best possible scenario, worst conditions
%Near extinction starting population - [E_0 = 0; M_0 = 1; F_0 = 1]
%Minimal R_0
%Low K (environmental carrying capacity)
%Low v (probability of emergence)
clear

%System parameters
K = 5000;
B_E = 7.46;
v_E = 0.005;
d_E = 0.046;
d_F = 0.046;
d_M = 0.139;
v = 0.45;
eta = 0.7;

%Spacing, day amount, f(x)
h = 0.001;
D = 400;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * ((eta * x(2)) / (1 + eta * x(2))) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3)];

%Solution and R_0
x = RK4(f, h, D, [0; 1; 1]);
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

%Plot of E, M and F
T = 0:h:D;
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
title('Near extinction starting population, worst conditions');
pbaspect([2 1 1]); 
hold off;

%Initial condition: worst possible scenario, near perfect conditions
%High starting population - [E_0 = 10000; M_0 = 2000; F_0 = 5000]
%Maximal R_0
%High K (environmental carrying capacity)
%High v (probability of emergence)
clear

%System parameters
K = 30000;
B_E = 14.85;
v_E = 0.25;
d_E = 0.023;
d_F = 0.033;
d_M = 0.077;
v = 0.55;
eta = 0.7;

%Spacing, day amount, f(x)
h = 0.001;
D = 400;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * ((eta * x(2)) / (1 + eta * x(2))) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3)];

%Solution and R_0
x = RK4(f, h, D, [10000; 2000; 5000]);
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

%Plot of E, M and F
T = 0:h:D;
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
title('High starting population, near perfect conditions');
pbaspect([2 1 1]); 
hold off;

%Initial condition: recent invasion/establishment, likely conditions
%Modest starting population - [E_0 = 50; M_0 = 10; F_0 = 20]
%Moderate R_0
%Moderate K (environmental carrying capacity)
%Likely v (probability of emergence)
clear

%System parameters
K = 22200;
B_E = 10;
v_E = 0.05;
d_E = 0.03;
d_F = 0.04;
d_M = 0.1;
v = 0.49;
eta = 0.7;

%Spacing, day amount, f(x)
h = 0.001;
D = 400;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * ((eta * x(2)) / (1 + eta * x(2))) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3)];

%Solution and R_0
x = RK4(f, h, D, [50; 10; 20]);
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

%Plot of E, M and F
T = 0:h:D;
figure;
plot(T, x(1,:), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T, x(2,:), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T, x(3,:), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
%xlim([0 20])
%ylim([0 500]), grid on
xlabel('Time (days)');
ylabel('Population');
legend;
title('Recent invasion/establishment, likely conditions');
pbaspect([2 1 1]); 
hold off;

%Initial condition: likely conditions
%Near extinction starting population - [E_0 = 0; M_0 = 1; F_0 = 1]
%Moderate R_0
%Moderate K (environmental carrying capacity)
%Likely v (probability of emergence)
clear

%System parameters
K = 22200;
B_E = 10;
v_E = 0.05;
d_E = 0.03;
d_F = 0.04;
d_M = 0.1;
v = 0.49;
eta = 0.7;

%Spacing, day amount, f(x)
h = 0.001;
D = 400;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * ((eta * x(2)) / (1 + eta * x(2))) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3)];

%Solution and R_0
x = RK4(f, h, D, [0; 1; 1]);
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

%Plot of E, M and F
T = 0:h:D;
figure(4);
plot(T, x(1,:), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T, x(2,:), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T, x(3,:), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
%xlim([0 10])
%ylim([0 500]), grid on

xlabel('Time (days)');
ylabel('Population');
legend;
title('Near extinction starting population, likely conditions');
pbaspect([2 1 1]); 
hold off;