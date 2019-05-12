function [dist, cone_dist] = get_zeroCondeDistanceParametrization(sz, d)

    [x, y, z] = get_basisGrid(sz);

    B = d(1) * x + d(2) * y + d(3) * z;
    R = sqrt(x.^2 + y.^2 + z.^2);
    cosTheta = abs(B ./ R);
    rad2deg = 180 / pi;
    theta_rad = acos(cosTheta);

    angles_deg = 0:90;
    magic_angle_rad = acos(1/sqrt(3));

    cone1 = abs(round(R .* sin(magic_angle_rad + theta_rad)));
    cone2 = abs(round(R .* sin(magic_angle_rad - theta_rad)));
    conemask = zeros(sz);
    conemask(cone1 < cone2) = 1;
    cone_dist = cone1 .* conemask + cone2 .* ~conemask;

    dist = 0:max(cone_dist(:));

end