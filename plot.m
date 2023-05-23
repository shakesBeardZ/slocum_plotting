clear
close all

addpath('gsw_matlab')
addpath('export_fig')

Zcbar = load('Zaro_colorbars.mat');
vars  = {'T'; 'C'; 'S'; 'chla'; 'CDOM'; 'BB'; 'DO'};
varN  = {'Temperature'; 'conductivity'; 'Salinity'; 'chlorophyl a'; 'CDOM'; 'Backscatter'; 'Dissolved Oxygen'};
units = {'^oC', 'S m^-1', 'PSU', };
clims = [[21 29]; [4 7]; [39 40.5]; [0 1]; [0.2 .7]; [0 3e-4]; [20 180]];
cbars = {cmocean('thermal'), Zcbar.ctemp2, cmocean('haline'), cmocean('algae'), cmocean('matter'), cmocean('matter'), 'parula'};

outfolder = './png';
ncfile = '/home/ubuntu/gliders/dbdreader/processed/unit_907-2023-133-4_timeseries.nc';


Time  = ncread(ncfile, 'time') / (3600 * 24) + datenum(1970, 1, 1);
T     = ncread(ncfile, 'temperature');
C     = ncread(ncfile, 'conductivity');
S     = ncread(ncfile, 'salinity');
P     = ncread(ncfile, 'press');
depth = ncread(ncfile, 'depth');
chla  = ncread(ncfile, 'chla');
DO    = ncread(ncfile, 'DO');
CDOM  = ncread(ncfile, 'CDOM');
BB    = ncread(ncfile, 'BB');
Dens  = ncread(ncfile, 'dens');
Lat   = ncread(ncfile, 'lat');
Lon   = ncread(ncfile, 'lon');

%%
figure
yyaxis('left')
plot(Time, T)
yyaxis('right')
plot(Time,depth)
ylim([-500 0]);
datetick('x','dd/mm', 'keeplimits');
%%
Timestep = (max(Time) - min(Time)) / 1000;
TimeBins = min(Time):Timestep:max(Time);
Depths = -1:-1:-600;

[XX, YY] = meshgrid(TimeBins, Depths);

ind = find(depth > 0 | isnan(depth) | depth < -1000);
depth(ind) = [];
Time(ind)  = [];
T(ind)     = [];
C(ind)     = [];
S(ind)     = [];
P(ind)     = [];
chla(ind)  = [];
DO(ind)    = [];
CDOM(ind)  = [];
BB (ind)   = [];
Dens(ind)  = [];
Lat(ind)   = [];
Lon(ind)   = [];

F = scatteredInterpolant(Time, depth, Dens, 'natural', 'none');
Dbin = F(XX, YY);

%%
TSplot(T, S, depth, [outfolder '/T-S'])

%%
for i=1:length(vars)
    var = vars{i};
    F = scatteredInterpolant(Time, depth, eval(var), 'natural', 'none');
    Vbin = F(XX, YY);
    plotHovm(XX, YY, Vbin, Dbin, clims(i,:), varN{i}, varN{i}, cbars{i}, [outfolder '/' varN{i}])
end

%%
function TSplot(T, S, depth, filename)
    figure('Position', [50 100 800 600])
    set(gcf,'color','w');
    set(gca,'fontsize', 18);
    hold on; box on

    scatter(S, T, 1, -depth, "." )
    cb = colorbar;
    cb.Label.String = "Depth";
    caxis([0 600])
    SS = 38.6:0.001:41;
    TT = 21:0.01:29;
    [Tp, Sp] = meshgrid(TT,SS);

    PP = zeros(size(Tp));
    D = gsw_rho(Sp,Tp,PP);
    contour(Sp,Tp,D, 'EdgeColor','#666666', 'ShowText',"on")
    xlabel("Salinity (PSU)")
    ylabel("Temeprature (^oC)")
    xlim([38.5 41])
    if strcmp(filename, '') == 0
    % savefig(filename)
        export_fig([filename '.png'], '-m2')
    end
end


%%
function plotHovm(lm, tm, Vm, Dens, climits, cLabel, Title, cbar, filename)
DensLevels = [1025 1026 1027 1028 1028.5];

figure('Position', [50 100 800 600])
set(gcf,'color','w');
set(gca,'fontsize', 18);
hold on; box on

h = pcolor(lm, tm, Vm);
set(h, 'edgecolor','none')
colormap(cbar);
shading interp

datetick('x','dd/mm', 'keeplimits');

ylim([-600 0])
if ~isnan(climits(1,1))
    caxis(climits);
end
cb = colorbar;
cb.Label.String = cLabel;

[C,h] = contour(lm, tm, Dens, DensLevels, 'LineWidth', 2, 'EdgeColor','w');
clabel(C,h,'FontSize',10,'Color','w')

title(Title);

if strcmp(filename, '') == 0
    % savefig(filename)
    export_fig([filename '.png'], '-m2')
end

end
