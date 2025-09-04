%SIT model with feedback law depending only on aquatic phase mosquitoes%

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

%Determining omega (w) lower bound
R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));

%r = 1 + (2 * d_M / (eta * K * (1 - v) * v_E)) * (1 + sqrt(1 + eta * K * (1 - v) * v_E / d_M));

w_lower_bound = (R_0 - 1) * (d_s * (1 - v) * v_E) / (d_M * y_s * 1); %previously 2.3122, now 2.2390

w = 4.5;

%Spacing, day amount, f(x)
h = 0.001;
D = 800;
f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            w * x(1) - d_s * x(4)];

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
u_plot = w * x(1, :);
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
plot(T(1:25000), MM_plot(1:25000), 'Color', [0.6350 0.0780 0.1840], 'DisplayName', ...
    'Michaelis-Menten factor', LineWidth=1.5);
xlabel('Time (days)');
ylabel('Michaelis-Menten factor');
legend;
title('Michaelis-Menten factor');
pbaspect([2 1 1]); 

%Sterilized male mosquitoes released in the 1st 800 days (aquatic phase
%feedback)
TM_s = trapz(T, u_plot);
fprintf('Sterilized male mosquitoes released in the 1st 800 days (aquatic phase feedback) = %d\n', TM_s)
%%
%Finding w that guarantees extiction and minimizes total sterile males

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
w = 2.5;

%Storing integral of u for S simulations
S = 10;
res = zeros(S, 3);

for i = 1:S
    %Spacing, day amount, f(x)
    h = 0.001;
    D = 800;
    T = 0:h:D;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
                (1-v) * v_E * x(1) - d_M * x(2);
                v * v_E * x(1) - d_F * x(3);
                w * x(1) - d_s * x(4)];
        
    %Solution
    x = RK4_mod(f, h, 1000, [21910; 5587; 13419; 0]);
    D = (length(x)-1) * h;
    T = 0:h:D;
    res(i,3) = trapz(T, w * x(1, :));
    res(i,2) = D;
    res(i,1) = w;

    w = w + 0.5;
end

%Plot of TM_s in function of lambda
xx = 2.5:.25:w;
yy = spline(res(:, 1),res(:, 3),xx);
plot(res(:, 1), res(:, 3), 'o', 'MarkerEdgeColor', [0 0 0.5], 'LineWidth', 1, 'DisplayName', '\omega');
hold on
plot(xx, yy, 'Color', [0.3010 0.7450 0.9330], ...
     'DisplayName', 'Total M_s released', 'LineWidth', 1.2);
xlabel('\omega');
ylabel('Total M_s Released');
legend;
title('Evolution of the control cost');
pbaspect([2 1 1]); 
xlim([2, 7])
hold off
%%
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
w = 4;

%Storing integral of u for S simulations
S = 10;
TM_s_diff_w = zeros(S, 1);
EMFSum = zeros(S, 1);

for i = 1:S
    %Spacing, day amount, f(x)
    h = 0.001;
    D = 800;
    T = 0:h:D;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
                (1-v) * v_E * x(1) - d_M * x(2);
                v * v_E * x(1) - d_F * x(3);
                w * x(1) - d_s * x(4)];
    
    %Solution
    x = RK4(f, h, D, [21910; 5587; 13419; 0]);

    TM_s_diff_w(i) = trapz(T, w * x(1, :));
    EMFSum(i) = sum(x(1:3, 800001));

    w = w + 0.1;
end
%%
%Robustness test (normal não ser robusto porque limite inferior de omega
%para que haja extinção depende dos parâmetros)
clear

h = 0.001;
D = 900;
S = 200; %should be 200
res = zeros(S, 1 + D / h);
%
K = 22200;
v = 0.49;
w = 4.3;
eta = 0.7;

w_LB = zeros(S, 1);
B_E_values = zeros(S, 1);
v_E_values = zeros(S, 1);
d_E_values = zeros(S, 1);
d_F_values = zeros(S, 1);
d_M_values = zeros(S, 1);
d_s_values = zeros(S, 1);
y_s_values = zeros(S, 1);

tic;
for i = 1:S
    B_E = unifrnd(7.46, 14.85);
    v_E = unifrnd(0.005, 0.25);
    d_E = unifrnd(0.023, 0.046);
    d_F = unifrnd(0.033, 0.046);
    d_M = unifrnd(0.077, 0.139);
    d_s = unifrnd(0.077, 0.139);
    y_s = unifrnd(0.5, 1);

    B_E_values(i) = B_E;
    v_E_values(i) = v_E;
    d_E_values(i) = d_E;
    d_F_values(i) = d_F;
    d_M_values(i) = d_M;
    d_s_values(i) = d_s;
    y_s_values(i) = y_s;

    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            w * x(1) - d_s * x(4)];

   %x = RK4(f, h, D, [normrnd(K, (0.1 * K) / 3);
   % normrnd(5587, (0.1 * 5587) / 3); normrnd(13419, (0.1 * 13419) / 3); 0]);
    x = RK4(f, h, D, [normrnd(K, (0.8 * K) / 3); normrnd(K, (0.8 * K) / 3); normrnd(K, (0.8 * K) / 3); normrnd(K, (0.8 * K) / 3)]);
    res(i, :) = sum(x(1:3, :));

    R_0 = (B_E * v * v_E) / (d_F * (v_E + d_E));
    r = 1 + (2 * d_M / (eta * K * (1 - v) * v_E)) * (1 + sqrt(1 + eta * K * (1 - v) * v_E / d_M));

    w_LB(i) = (R_0 - r) * (d_s * (1 - v) * v_E) / (d_M * y_s * r);
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

% figure;
% plot(1:1:20, B_E_values)
% xlabel('Simulation')
% ylabel('B_E')
% 
% figure;
% plot(1:1:20, v_E_values)
% xlabel('Simulation')
% ylabel('v_E')
% 
% figure;
% plot(1:1:20, d_E_values)
% xlabel('Simulation')
% ylabel('d_E')
% 
% figure;
% plot(1:1:20, d_F_values)
% xlabel('Simulation')
% ylabel('d_F')
% 
% figure;
% plot(1:1:20, d_M_values)
% xlabel('Simulation')
% ylabel('d_M')
% 
% figure;
% plot(1:1:20, d_s_values)
% xlabel('Simulation')
% ylabel('d_s')
% 
% figure;
% plot(1:1:20, y_s_values)
% xlabel('Simulation')
% ylabel('y_s')
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
eta = 0.7;
%w = 2.23;

%Storing for S simulations
S = 10;
res = zeros(S, 3);

% %First simulation
%     h = 0.001;
%     f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
%             (1-v) * v_E * x(1) - d_M * x(2);
%             v * v_E * x(1) - d_F * x(3);
%             w * x(1) - d_s * x(4)];
% 
%     %Solution
%     x = RK4_mod(f, h, 6000, [21910; 5587; 13419; 0]);
%     D = (length(x)-1) * h;
%     T = 0:h:D;
%     res(1,3) = trapz(T, w * x(1, :));
%     res(1,2) = D;
%     res(1,1) = w;

%Remaining simulations
w = 2.5;
for i = 1:S
    %Spacing, day amount, f(x)
    h = 0.001;
    f = @(t, x) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
            (1-v) * v_E * x(1) - d_M * x(2);
            v * v_E * x(1) - d_F * x(3);
            w * x(1) - d_s * x(4)];
    
    %Solution
    x = RK4_mod(f, h, 1000, [21910; 5587; 13419; 0]);
    D = (length(x)-1) * h;
    T = 0:h:D;
    res(i,3) = trapz(T, w * x(1, :));
    res(i,2) = D;
    res(i,1) = w;

    w = w + 0.5;
end