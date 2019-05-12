function metrics = compute_1dMetrics(img1, img2, options)

    if nargin < 3
        options = struct('direction', [0, 0, 1], ...
                         'smoothing', 0.05);
    end
    metrics.options = options;

    options.parametrization = 'shell';
    [metrics.x_shell, metrics.y_shell] = esp(img1, img2, options);

    options.parametrization = 'cone';
    [metrics.x_cone, metrics.y_cone] = esp(img1, img2, options);

    options.parametrization = 'zeroConeDistance';
    [metrics.x_zcd, metrics.y_zcd] = esp(img1, img2, options);

end