function result = Newton_Raphson(f, initial_guess)
    tolerance = 1E-6;
    x0 = initial_guess;
    error = 10000;

    while error > tolerance
        d_f = (f(x0 + tolerance) - f(x0 - tolerance)) / (2 * tolerance);
        x_n = x0 - f(x0) / d_f;
        error = abs(x_n - x0);
        x0 = x_n;
    end

    result = x0;
end