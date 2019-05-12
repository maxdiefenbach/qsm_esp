function h = plot_metrics(metricsList)

    h = figure;

    legend_str = {};
    for i  = 1:length(metricsList)
            metrics = metricsList(i)
            legend_str{i} = show_scalarMetrics(metrics);

            subplot(1, 3, 1);
            hold on
            plot(metrics.x_shell, metrics.y_shell);
            axis tight;
            xlabel('radius');

            subplot(1, 3, 2);
            hold on;
            plot(metrics.x_cone, metrics.y_cone);
            axis tight;
            xlabel('azimuth');

            subplot(1, 3, 3);
            hold on
            plot(metrics.x_zcd, metrics.y_zcd);
            axis tight;
            xlabel('zero-cone distance');

        end

    subplot(1, 3, 1);
    axis tight;
    xlabel('radius');

    subplot(1, 3, 2);
    YLIM = get(gca, 'YLIM');
    plot([54.7, 54.7], [0, YLIM(2)], 'black');

    axis tight;
    xlabel('aximuth');

    subplot(1, 3, 3);
    axis tight;
    xlabel('zero-cone distance');
    legend(legend_str)

    % annotation('textbox', [0.5, .9, .1, .1], ...
    %            'String', ...
    %            sprintf('rmse = %f, hfen = %f, ssim = %f', ...
    %                    metrics.rmse, metrics.hfen, metrics.ssim), ...
    %            'EdgeColor','none');

end