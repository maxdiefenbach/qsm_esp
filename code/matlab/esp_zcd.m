%% Description: function to compute zero-cone-distance-esp
%%
%% Diefenbach M. N., Boehm C., Meineke J., Liu C., Karampions D. C.
%% "One-Dimensional k-Space Metrics on Cone Surfaces for Quantitative Susceptibility Mapping",
%% Proceedings 27. Annual Meeting International Societz for Magnetic Resonance in Medicine, Montreal, 2019
%% Oral presentation: Monday, 13 May 2019, QSM & ETM Session, 4pm - 6pm
%% Abstract ID: 3839
%%
%% Authors: Christof Boehm, Maximillian N. Diefenbach
%%
%% Date created: February 18, 2019

function [x, f, y] = esp_zcd(refF2, errormapF2, d)
    FOV = size(refF2);

    [xm, ym, zm] = get_grid(FOV);

    R = sqrt(xm.^2 + ym.^2 + zm.^2);
    B = d(1) * xm + d(2) * ym + d(3) * zm;

    Theta = acos(B ./ R);
    magic_angle_rad = acos(1/sqrt(3));

    cone1 = abs(round(R .* sin(magic_angle_rad + Theta)));
    cone2 = abs(round(R .* sin(magic_angle_rad - Theta)));
    conemask = zeros(FOV);
    conemask(cone1 < cone2) = 1;
    cone = cone1 .* conemask + cone2 .* ~conemask;

    max_dist = max(cone(:));
    f = zeros(max_dist + 1, 1);
    y = zeros(max_dist + 1, 1);

    for i = 0:max_dist
        mask = cone == i;
        f(i+1) = sum(errormapF2(mask));
        y(i+1) = sum(refF2(mask));
    end

    x = 0:max_dist;
end