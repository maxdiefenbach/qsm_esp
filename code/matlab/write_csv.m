function write_csv(hfig, nametag)

    it = 0;

    fcs = get(hfig, 'Children');
    for i = numel(fcs):-1:1

        if strcmp(class(fcs(i)), 'matlab.graphics.axis.Axes')

            axis = fcs(i);
            acs = get(axis, 'Children'); % sorted backwards

            for j = numel(acs):-1:1

                line = acs(j);
                XData = get(line, 'XData');
                YData = get(line, 'YData');

                t = table(XData(:), YData(:));
                it = it + 1;
                writetable(t, [nametag '_table' num2str(it) '.csv']);

% t.(['xVar' num2str(it)]) = XData(:);
% t.(['yVar' num2str(it)]) = YData(:);

            end

        end

    end
end