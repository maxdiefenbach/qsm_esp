%% Description: computing line-plots of Figure 3 of the corresponding conference abstract
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

%% Load data
clear; clc;
load('../../data/reconchallenge2016/data/chi_l2.mat')
chi{1} = chi_l2;
clear chi_l2
string_plot{1} = 'L2';

load('../../data/reconchallenge2016/data/chi_tkd.mat')
chi{2} = chi_tkd;
string_plot{2} = 'TKD';
clear chi_tkd

load('../../data/reconchallenge2016/data/chi_cosmos.mat')
chi{3} = chi_cosmos;
string_plot{3} = 'COSMOS';
clear chi_cosmos

%% load STI data as reference
load('../../data/reconchallenge2016/data/chi_33.mat')

%% run all ESPs
B0dir = [0, 0, 1];
for i = 1:3
    [x_shell{i}, y_shell{i}] = esp_shell(chi_33, chi{i});
    [x_zcd{i}, y_zcd{i}] = esp_zcd(chi_33, chi{i}, B0dir);
    [x_cone{i}, y_cone{i}] = esp_cone(chi_33, chi{i}, B0dir);
end

%% Plot
figure;
hold on
for i = 1:3
    plot(x_shell{i}, y_shell{i}, 'LineWidth', 2);
end
xlabel('radii')
xlim([0 max(x_shell{1}(:))])
legend(string_plot, 'FontSize', 14)
hold off
grid on

figure;
hold on
for i = 1:3
    plot(x_cone{i}, y_cone{i}, 'LineWidth', 2);
end
xlabel('azimuth')
legend(string_plot, 'FontSize', 14)
hold off
grid on

figure;
hold on
for i = 1:3
    plot(x_zcd{i}, y_zcd{i}, 'LineWidth', 2);
end
xlabel('distance to zero cone surface')
xlim([0 max(x_zcd{1}(:))])
legend(string_plot, 'FontSize', 14, 'Location', 'northwest')
hold off
grid on