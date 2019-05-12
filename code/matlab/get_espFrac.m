function [x_esp, y_esp] = get_espFrac(x, f, y, smoothing)

    if nargin < 4
        smoothing = 1
    end

    x_esp = x;

    if smoothing >= 0
        spl_f = fit(x(:), f(:), 'smoothingspline', 'SmoothingParam', smoothing);
        spl_y = fit(x(:), y(:), 'smoothingspline', 'SmoothingParam', smoothing);
        y_esp = sqrt(spl_f(x) ./ spl_y(x));

    else
        y_esp = sqrt(f ./ y);

    end

    y_esp = y_esp(:).';

end