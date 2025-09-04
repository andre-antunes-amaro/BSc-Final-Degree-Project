function x = RK4_mod(f, h, t, VI)

    T = 0:h:t;
    x = zeros(length(VI), length(T));
    x(:,1) = VI;
    K = 22200;
    
    for i = 1:(length(T)-1)
        F1 = f(T(i), x(:, i));

        F2 = f(T(i) + h/2, x(:, i) + (h/2) .* F1);

        F3 = f(T(i) + h/2, x(:, i) + (h/2) .* F2);

        F4 = f(T(i) + h, x(:, i) + h .* F3);
        if x(1, i) > K / 100
            x(:,i+1) = x(:,i) + (h/6) .* (F1 + 2 * F2 + 2 * F3 + F4);
        else
            x = x(:, 1:i);
            break
        end
    end
end