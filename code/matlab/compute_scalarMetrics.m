function scalarMetrics = compute_scalarMetrics(img1, img2)

% img1 = true
    scalarMetrics.rmse = compute_rmse(img2, img1);
    scalarMetrics.hfen = compute_hfen(img2, img1);
    scalarMetrics.ssim = compute_ssim(img2, img1);

end