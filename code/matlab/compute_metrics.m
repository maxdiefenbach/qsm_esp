function metrics = compute_metrics(img1, img2, options)
% img1 = true

    if nargin < 3
        options = struct('direction', [0, 0, 1], ...
                         'smoothing', 0.05);
    end

    metrics = merge_Structs(compute_scalarMetrics(img1, img2), ...
                            compute_1dMetrics(img1, img2, options));

end