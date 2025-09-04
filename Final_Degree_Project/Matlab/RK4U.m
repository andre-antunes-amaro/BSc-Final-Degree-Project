function x = RK4U(f, u, h, D, VI)
    N = D/h + 1;

    x = zeros(length(VI), N);
    x(:,1) = VI;
    
    for i = 1:(N-1)
        F1 = f(x(:, i), u(i));

        F2 = f(x(:, i) + (h/2) .* F1, (u(i) + u(i+1))/2);

        F3 = f(x(:, i) + (h/2) .* F2, (u(i) + u(i+1))/2);

        F4 = f(x(:, i) + h .* F3, u(i+1));

        x(:,i+1) = x(:,i) + (h/6) .* (F1 + 2 * F2 + 2 * F3 + F4);
    end
end