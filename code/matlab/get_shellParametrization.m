function [radii, R] = get_shellParametrization(sz)

    [x, y, z] = get_basisGrid(sz);

    R = sqrt(x.^2 + y.^2 + z.^2);

    radii = 1:(max(sz) / 2 + 1);

end