function [xm, ym, zm] = get_grid(FOV)

    xbase = linspace(-FOV(1)/2, FOV(1)/2 - 1, FOV(1));
    ybase = linspace(-FOV(2)/2, FOV(2)/2 - 1, FOV(2));
    zbase = linspace(-FOV(3)/2, FOV(3)/2 - 1, FOV(3));
    [xm, ym, zm] = ndgrid(xbase, ybase, zbase);

end