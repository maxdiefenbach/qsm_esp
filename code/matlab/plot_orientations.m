function h = plot_orientations(array3d, varargin)

    center = ceil((size(array3d) + 1) / 2);
    ix = center(1);
    iy = center(2);
    iz = center(3);
    clim = [min(array3d(:)), max(array3d(:))];
    for v = 1:(nargin-1)
        v
        nargin
        if strcmp(varargin{v}, 'ix')
            ix = varargin{v+1};
        end
        if strcmp(varargin{v}, 'iy')
            iy = varargin{v+1};
        end
        if strcmp(varargin{v}, 'iz')
            iz = varargin{v+1};
        end
        if strcmp(varargin{v}, 'CLim')
            clim = varargin{v+1};
        end
    end
    arrK = squeeze(array3d(:, :, iz));
    arrJ = permute(squeeze(array3d(:, iy, :)), [2, 1, 3]);
    arrI = permute(squeeze(array3d(ix, :, :)), [2, 1, 3]);

    h = figure;
    colormap viridis;

    subplot(1, 3, 1)
    imagesc(arrK)
    hold on
    plot(iy, ix, '+', 'Color', 'red', 'markers', 20)

    axis off

    subplot(1, 3, 2)
    imagesc(arrJ)
    hold on
    plot(iz, ix, '+', 'Color', 'red', 'markers', 20)
    axis off

    subplot(1, 3, 3)
    imagesc(arrI)
    hold on
    plot(iz, iy, '+', 'Color', 'red', 'markers', 20)
    colorbar
    axis off

    children = get(h, 'Children')
    for i = 1:numel(children)
        if strcmp(class(children(i)), 'matlab.graphics.axis.Axes')
            set(children(i), 'CLim', clim)
        end
    end

end