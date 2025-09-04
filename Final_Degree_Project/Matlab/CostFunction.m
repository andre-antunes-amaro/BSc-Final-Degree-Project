function J = CostFunction(X)
    a1 = 1;
    a2 = 0;
    a3 = 0;
    B = 1;
    h = 10^-1;
    
    N = length(X)/5;

    E = X(1:N);
    M = X(N+1:2*N);
    F = X(2*N+1:3*N);
    u = X(4*N+1:5*N);

    % J = h/3 * ( (a1*F(1)^2 + B*u(1)) + ...
    %             4*sum(a1*F(2:2:end-1).^2 + B*u(2:2:end-1)) + ...
    %             2*sum(a1*F(3:2:end-2).^2 + B*u(3:2:end-2)) + ...
    %             (a1*F(end)^2 + B*u(end)) );

    J = h/3 * ( (a1*E(1)^2 + a2*M(1)^2 + a3*F(1)^2 + B*abs(u(1))) + ...
                4*sum(a1*E(2:2:end-1).^2 + a2*M(2:2:end-1).^2 +a3*F(2:2:end-1).^2 + B*abs(u(2:2:end-1))) + ...
                2*sum(a1*E(3:2:end-2).^2 + a2*M(3:2:end-2).^2 +a3*F(3:2:end-2).^2 + B*abs(u(3:2:end-2))) + ...
                (a1*E(end)^2 + a2*M(end)^2 +a3*F(end)^2 + B*abs(u(end))) );

    % J = sum((u - 10^5).^2);
end