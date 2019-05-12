function img12 = fuse_zcd(direction, img1, img2, metrics1, metrics2)


    FOV = size(img1);
    direction = [0, 0, 1];
    [x, y, z] = get_grid(FOV);
    R = sqrt(x.^2 + y.^2 + z.^2);
    B = direction(1) * x + direction(2) * y + direction(3) * z;
    Theta = acos(B ./ R);
    magic_angle_rad = acos(1/sqrt(3));
    cone1 = abs(round(R .* sin(magic_angle_rad + Theta)));
    cone2 = abs(round(R .* sin(magic_angle_rad - Theta)));
    conemask = zeros(FOV);
    conemask(cone1 < cone2) = 1;
    cone = cone1 .* conemask + cone2 .* ~conemask;

    take = metrics2.y_zcd <= metrics1.y_zcd;
    Img_fused = fft3c(img1);
    Img_maybe = fft3c(img2);
    for ir = 1:length(metrics1.x_zcd)
        zcd = metrics1.x_zcd(ir);
        mask_zcd = cone == zcd;
        if take(ir)
            Img_fused(mask_zcd) = Img_maybe(mask_zcd);
        end
    end
    img12 = ifft3c(Img_fused);

end