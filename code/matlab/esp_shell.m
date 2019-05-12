%% Description: function to compute shell-esp
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

function [x, f, y] = esp_shell(refF2, errormapF2)
    FOV = size(refF2);

    [xm, ym, zm] = get_grid(FOV);

    R = sqrt(xm.^2 + ym.^2 + zm.^2);

    radii = 1:(max(FOV) / 2) + 1;

    f = zeros(length(radii), 1);
    y = zeros(length(radii), 1);

    for i = radii
        mask = ((i-1) <= R) & (R < i);
        f(i) = sum(errormapF2(mask));
        y(i) = sum(refF2(mask));
    end

    x = radii;
end