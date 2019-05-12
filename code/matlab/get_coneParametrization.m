function [angles_deg, theta_deg] = get_coneParametrization(sz, d)

    [x, y, z] = get_basisGrid(sz);

    B = d(1) * x + d(2) * y + d(3) * z;
    R = sqrt(x.^2 + y.^2 + z.^2);
    cosTheta = abs(B ./ R);
    rad2deg = 180 / pi;

    theta_deg = round(rad2deg .* acos(cosTheta));
    angles_deg = 0:90;

end