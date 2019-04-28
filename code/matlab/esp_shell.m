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

function [x, y] = esp_shell(ref, img)
    FOV = size(ref);
    refF = ftn(ref);
    imgF = ftn(img);
    refF_squared = abs(refF) .^ 2;
    errormap = abs(imgF - refF) .^ 2;

    xbase = linspace(-FOV(1)/2, FOV(1)/2 - 1, FOV(1));
    ybase = linspace(-FOV(2)/2, FOV(2)/2 - 1, FOV(2));
    zbase = linspace(-FOV(3)/2, FOV(3)/2 - 1, FOV(3));
    [xm, ym, zm] = ndgrid(xbase, ybase, zbase);

    R = sqrt(xm .^ 2 + ym .^ 2 + zm .^ 2);

    radii = 1:(max(FOV) / 2) + 1;

    error1 = zeros(length(radii), 1);
    error2 = zeros(length(radii), 1);

    for i = radii
        mask = ((i-1) <= R) & (R < i);
        tmp_errormap = errormap .* mask;
        tmp_refF_squared = refF_squared .* mask;

        error1(i) = sum(tmp_errormap(:));
        error2(i) = sum(tmp_refF_squared(:));
    end

    x = radii;
    y = sqrt(error1 ./ error2);
    y = y(:).';
end