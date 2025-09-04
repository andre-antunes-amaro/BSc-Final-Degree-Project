function [c, ceq] = Constraints(X)
    h = 10^-1;
    D = 30;

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
    f = @(x, u) [B_E * x(3) * (1 - x(1) / K) * (eta * x(2) / (1 + eta * (x(2) + y_s * x(4))) ) - (v_E + d_E) * x(1);
                (1-v) * v_E * x(1) - d_M * x(2);
                v * v_E * x(1) - d_F * x(3);
                u - d_s * x(4)];

    N = D/h + 1;

    c = [];

    ceq = zeros(4*(N-1), 1);

    x = [X(1:N)'; X(N+1:2*N)'; X(2*N+1:3*N)'; X(3*N+1:4*N)'];
    x(:, 1) = [21910; 5587; 13419; 0];
    %u = zeros(1, N) + 10^5;
    u = X(4*N+1:5*N);
        
    j = 1;
    for i = 1:(N-1)

        F1 = f(x(:, i), u(i));

        F2 = f(x(:, i) + (h/2) .* F1, (u(i) + u(i+1))/2);

        F3 = f(x(:, i) + (h/2) .* F2, (u(i) + u(i+1))/2);

        F4 = f(x(:, i) + h .* F3, u(i+1));

        ceq(j:j+3) = x(:,i+1) - x(:,i) - (h/6) .* (F1 + 2 * F2 + 2 * F3 + F4);
        
        %ceq(j:j+3) = x(:,i+1) - x(:,i) - h.* f(x(:,i+1), u(i+1));

        j = j+4;
    end

    %c = x(1, end) - 10^3;
end

