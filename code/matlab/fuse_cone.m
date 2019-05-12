function img12 = fuse_cone(direction, img1, img2, metrics1, metrics2)

    FOV = size(img1);
    direction = [0, 0, 1];
    [x, y, z] = get_grid(FOV);
    R = sqrt(x.^2 + y.^2 + z.^2);
    B = direction(1) * x + direction(2) * y + direction(3) * z;
    Theta = acos(B ./ R);

    take = metrics2.y_cone <= metrics1.y_cone;
    Img_fused = fft3c(img1);
    Img_maybe = fft3c(img2);
    for ir = 1:length(metrics1.x_cone)
        angle_rad = metrics1.x_cone(ir)
        mask_cone = abs(B) == round(R .* sin(pi/2 - angle_rad));
        if take(ir)
            Img_fused(mask_cone) = Img_maybe(mask_cone);
        end
    end

    img12 = ifft3c(Img_fused);

end