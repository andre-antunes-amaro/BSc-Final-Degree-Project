function x = RK4(f, h, D, VI)

    T = 0:h:D;
    x = zeros(length(VI), length(T));
    x(:,1) = VI;
    
    for i = 1:(length(T)-1)
        F1 = f(T(i), x(:, i));

        F2 = f(T(i) + h/2, x(:, i) + (h/2) .* F1);

        F3 = f(T(i) + h/2, x(:, i) + (h/2) .* F2);

        F4 = f(T(i) + h, x(:, i) + h .* F3);

        x(:,i+1) = x(:,i) + (h/6) .* (F1 + 2 * F2 + 2 * F3 + F4);
    end
end