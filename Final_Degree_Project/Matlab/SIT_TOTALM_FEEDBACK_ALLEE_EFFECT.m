%SIT model (considering allee effect) with feedback law depending only on the total number of males%

%Initial condition: likely scenario for application 
%Starting population at persistance equilibrium -
%[E_0 = 21910; M_0 = 5587; F_0 = 13419; M_s_0 = 0]
%Moderate R_0
%Moderate K (environmental carrying capacity)
%Likely v (probability of emergence)
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

%k lower bound, considering R_1 (k) < r
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

r = 1 + (2 * d_M / (eta * K * (1 - v) * v_E)) * (1 + sqrt(1 + eta * K * (1 - v) * v_E / d_M));


k_lower_bound = d_s * (R_0 - r) / (R_0 + r * (y_s - 1)); %previously 0.11843, now 0.11838

k = 0.119;

%Spacing, day amount, f(x)
h = 0.001;
D = 4000;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            k * (x(2) + x(4)) - d_s * x(4)];

%Solution
x = RK4(f, h, D, [21910; 5587; 13419; 0]);

T = 0:h:D;
%Plot of E, M and F
figure;
plot(T(1:700/h), x(1,1:700/h), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T(1:700/h), x(2,1:700/h), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T(1:700/h), x(3,1:700/h), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of E, M and F');
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
u_plot = k * (x(2, :) + x(4, :));
figure;
plot(T, u_plot, 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Control function u', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Release rate');
legend;
title('Plot of the control u');
pbaspect([2 1 1]); 

%Plot of the MM factor
MM_plot = eta .* x(2, :) ./ (1 + eta .* (x(2, :) + y_s * x(4, :)));
figure;
plot(T(1:500000), MM_plot(1:500000), 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Michaelis-Menten factor', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Michaelis-Menten factor');
legend;
title('Michaelis-Menten factor');
pbaspect([2 1 1]); 

%Sterilized male mosquitoes released in the 1st 700 days (total males
%feedback)
TM_s = trapz(T(1:700/h), u_plot(1:700/h));
fprintf('Sterilized male mosquitoes released in the 1st 700 days (total males feedback) = %d\n', TM_s)
%%
% %Finding k that guarantees extiction and minimizes total sterile males
% 
% clear
% 
% %System parameters
% K = 22200;
% B_E = 10;
% y_s = 1;
% v_E = 0.05;
% d_E = 0.03;
% d_F = 0.04;
% d_M = 0.1;
% d_s = 0.12;
% v = 0.49;
% eta = 0.7;
% k = 0.1175;
% 
% %Storing integral of u for S simulations
% S = 14;
% TM_s_diff_k = zeros(S, 1);
% EMFSum = zeros(S, 1);
% 
% for i = 1:S
%     %Spacing, day amount, f(x)
%     h = 0.001;
%     D = 700;
%     T = 0:h:D;
%     f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
%             (1-v) * v_E * x(1) - d_M * x(2);
%             v * v_E * x(1) - d_F * x(3);
%             k * (x(2) + x(4)) - d_s * x(4)];
% 
%     %Solution
%     x = RK4(f, h, D, [21910; 5587; 13419; 0]);
% 
%     TM_s_diff_k(i) = trapz(T, k * (x(2, :) + x(4, :)));
%     EMFSum(i) = sum(x(1:3, 700001));
% 
%     k = k + 0.0001;
% end
%%
%Fixed parameters, random initial conditions (Robustness test 1)
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
k = 0.119;
eta = 0.7;

%Spacing, day amount, f(x), simulation amount
h = 0.001;
D = 1000;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            k * (x(2) + x(4)) - d_s * x(4)];

S = 20;
res = zeros(S, 1 + D / h);

tic;
for i = 1:S

    %x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3); normrnd(5587, (0.1 * 5587) / 3); normrnd(13419, (0.1 * 13419) / 3); 0]);
    x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3)]);
    res(i, :) = sum(x(1:3, :));
end
runtime = toc;

%Plot
T = 0:h:D;

figure;
plot(T, res(1, :));
hold on;
for i = 2:S
    plot(T, res(i, :))
end
xlabel('Time (days)');
ylabel('E+F+M');
pbaspect([2 1 1]); 
hold off;
%%
%Robustness test 2
clear

%System parameters
K = 22200;
B_E = 11;
y_s = 1;
v_E = 0.08;
d_E = 0.046;
d_F = 0.033;
d_M = 0.12;
d_s = 0.139;
v = 0.49;
k = 0.119;
eta = 0.7;

%Spacing, day amount, f(x)
h = 0.001;
D = 800;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            k * (x(2) + x(4)) - d_s * x(4)];

%Solution
x = RK4(f, h, D, [21910; 5587; 13419; 0]);

T = 0:h:D;
%Plot of E, M and F
figure;
plot(T, x(1, :), 'Color', [0 0.4470 0.7410], 'DisplayName', ...
    'Mosquito eggs E', LineWidth=1.5);
hold on;
plot(T, x(2, :), 'Color', [0.8500 0.3250 0.0980], 'DisplayName', ...
    'Male mosquitoes M', LineWidth=1.5);
plot(T, x(3, :), 'Color', [0.4940 0.1840 0.5560], 'DisplayName', ...
    'Female mosquitoes F', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Population');
legend;
title('Plot of E, M and F');
pbaspect([2 1 1]); 
hold off;

R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

r = 1 + (2 * d_M / (eta * K * (1 - v) * v_E)) * (1 + sqrt(1 + eta * K * (1 - v) * v_E / d_M));


k_lower_bound = d_s * (R_0 - r) / (R_0 + r * (y_s - 1));