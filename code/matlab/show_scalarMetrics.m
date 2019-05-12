function s = show_scalarMetrics(metrics)

    s = sprintf('rmse\thfen\tssim\n%f\t%f\t%f', metrics.rmse, metrics.hfen, metrics.ssim);

end