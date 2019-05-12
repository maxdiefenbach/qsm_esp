function [f, y] = get_kSpaceProfiles(refF2, errormapF2, x, P)

    N = length(x);
    f = zeros(N, 1);
    y = zeros(N, 1);

    for i = 1:N
        mask = P == x(i);
        f(i) = sum(errormapF2(mask));
        y(i) = sum(refF2(mask));
    end

end