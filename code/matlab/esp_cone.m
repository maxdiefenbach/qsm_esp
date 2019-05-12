%% Description: function to compute cone-esp
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

function [x, f, y] = esp_cone(refF2, errormapF2, d)
    FOV = size(refF2);

    [xm, ym, zm] = get_grid(FOV);

    R = sqrt(xm .^ 2 + ym .^ 2 + zm .^ 2);
    B = d(1) * xm + d(2) * ym + d(3) * zm;

    angles = 0:90;
    f = zeros(length(angles), 1);
    y = zeros(length(angles), 1);
    deg2rad = pi / 180;

    for i = angles
        angle_rad = i * deg2rad;
        mask = abs(B) == round(R .* sin(pi/2 - angle_rad));
        f(i+1) = sum(errormapF2(mask));
        y(i+1) = sum(refF2(mask));
    end

    x = angles;
%     y = sqrt(f ./ y);
%     y = y(:).';
end