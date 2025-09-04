%SIT model with feedback law depending only on wild males%

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
lambda = 22;

%Spacing, day amount, f(x)
h = 0.001;
D = 400;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) * (x(2) / (x(2) + y_s * x(4))) - d_F * x(3);
            lambda * x(2) - d_s * x(4)];

%Solution
x = RK4(f, h, D, [21910; 5587; 13419; 0]);

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
u_plot = lambda * x(2, :);
figure;
plot(T, u_plot, 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Control function u', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Release rate');
legend;
title('Plot of the control u');
pbaspect([2 1 1]); 

%Sterilized male mosquitoes released in the 1st 400 days (wild males
%feedback)
TM_s = trapz(T, u_plot);
fprintf('Sterilized male mosquitoes released in the 1st 400 days (wild males feedback) = %d\n', TM_s)
%%
%SIT model with feedback law depending only on wild males (lambda = 11)%

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
lambda = 9.06;

%Spacing, day amount, f(x)
h = 0.001;
D = 2000;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) * (x(2) / (x(2) + y_s * x(4))) - d_F * x(3);
            lambda * x(2) - d_s * x(4)];

%Solution
x = RK4(f, h, D, [21910; 5587; 13419; 0]);

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
u_plot = lambda * x(2, :);
figure;
plot(T, u_plot, 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Control function u', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Release rate');
legend;
title('Plot of the control u');
pbaspect([2 1 1]); 
%%
%Comparison table

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
lambda = 9.06;

%Storing for S simulations
S = 14;
res = zeros(S, 3);

%First simulation
    h = 0.001;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) * (x(2) / (x(2) + y_s * x(4))) - d_F * x(3);
            lambda * x(2) - d_s * x(4)];
    
    %Solution
    x = RK4_mod(f, h, 5000, [21910; 5587; 13419; 0]);
    D = (length(x)-1) * h;
    T = 0:h:D;
    res(1,3) = trapz(T, lambda * x(2, :));
    res(1,2) = D;
    res(1,1) = lambda;

%Remaining simulations
lambda = 10;
for i = 2:S
    %Spacing, day amount, f(x)
    h = 0.001;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) * (x(2) / (x(2) + y_s * x(4))) - d_F * x(3);
            lambda * x(2) - d_s * x(4)];
    
    %Solution
    x = RK4_mod(f, h, 900, [21910; 5587; 13419; 0]);
    D = (length(x)-1) * h;
    T = 0:h:D;
    res(i,3) = trapz(T, lambda * x(2, :));
    res(i,2) = D;
    res(i,1) = lambda;

    lambda = lambda + 1;
end
%%
% %Fixed parameters, random initial conditions
% K = 22200;
% B_E = 10;
% y_s = 1;
% v_E = 0.05;
% d_E = 0.03;
% d_F = 0.04;
% d_M = 0.1;
% d_s = 0.12;
% v = 0.49;
% lambda = 22;
% 
% D = 1000;
% S = 20; 
% res = zeros(S, 1 + D / h);
% 
% tic;
% for i = 1:S
% 
%     %x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3); normrnd(5587, (0.1 * 5587) / 3); normrnd(13419, (0.1 * 13419) / 3); 0]);
%     %x = RK4(f, h, D, [unifrnd(0, 50000); unifrnd(0, 50000); unifrnd(0, 50000); 0]);
%     x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3)]);
%     res(i, :) = sum(x);
% end
% runtime = toc;
% 
% %Plot
% T = 0:h:D;
% 
% figure;
% plot(T, res(1, :));
% hold on;
% for i = 2:S
%     plot(T, res(i, :))
% end
% xlabel('Time (days)');
% ylabel('E+F+M');
% pbaspect([2 1 1]); 
% hold off;
%%
%Robustness test 
clear

h = 0.001;
D = 900;
S = 200; %should be 200
res = zeros(S, 1 + D / h);
%
K = 22200;
B_E = 10;
y_s = 1;
v_E = 0.05;
d_E = 0.03;
d_F = 0.04;
d_M = 0.1;
d_s = 0.12;
v = 0.49;
lambda = 22;

tic;
for i = 1:S
    % B_E = unifrnd(7.46, 14.85);
    % v_E = unifrnd(0.005, 0.25);
    % d_E = unifrnd(0.023, 0.046);
    % d_F = unifrnd(0.033, 0.046);
    % d_M = unifrnd(0.077, 0.139);
    % d_s = unifrnd(0.077, 0.139);
    % y_s = unifrnd(0.5, 1);

    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) * (x(2) / (x(2) + y_s * x(4))) - d_F * x(3);
            lambda * x(2) - d_s * x(4)];

   %x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3);
   % normrnd(5587, (0.1 * 5587) / 3); normrnd(13419, (0.1 * 13419) / 3); 0]);
   x = RK4(f, h, D, [unifrnd(0, 10 * K); unifrnd(0, 10 * K); unifrnd(0, 10 * K); unifrnd(0, 10 * K);]);
   %x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3); normrnd(K, (0.1 * K) / 3)]);
   res(i, :) = sum(x(1:3, :));
end
runtime = toc;

%Plot
T = 0:h:D;

figure;
plot(T, log(res(1, :)));
hold on;
for i = 2:S
    plot(T, log(res(i, :)))
end
plot(T, zeros(1, length(T)), LineWidth=2, Color = 'r')
xlabel('Time (days)');
ylabel('log(E+F+M)');
pbaspect([2 1 1]); 
hold off;



