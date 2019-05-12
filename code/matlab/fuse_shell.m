function img12 = fuse_shell(img1, img2, metrics1, metrics2)

    [x, y, z] = get_grid(size(img1));
    R = sqrt(x.^2 + y.^2 + z.^2);
    take = metrics2.y_shell <= metrics1.y_shell;
    Img_fused = fft3c(img1);
    Img_maybe = fft3c(img2);
    for ir = 1:length(metrics1.x_shell)
        r = metrics1.x_shell(ir);
        shell = ((r-1) <= R) & (R < r);
        if take(ir)
            Img_fused(shell) = Img_maybe(shell);
        end
    end
    img12 = ifft3c(Img_fused);

end