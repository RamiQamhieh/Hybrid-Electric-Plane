% first run aircraftDatabase.m to load in data
aircraftDatabase

% plotting parameters
m_length = 500;
tstep = 30/3600;%15/3600; %[h] 15 sec -> hrs; time step size


%% Our aircraft
c = Aircraft;
c.Name = "Proposed Hybrid";
c.tstep = tstep
c.mission_(m_length)

%% Our aircraft - figures
figure;%('DefaultAxesFontSize',14);
% plot(c.r_time_elapsed,c.r_M)
hold on;
plot(c.r_time_elapsed,c.r_Mfuel)
ylabel('Fuel Mass [kg]');
xlabel('Time [h]');
% yline(c.)
h = min(c.r_Mfuel)-100;
xline(c.r_time_elapsed(find(c.r_stage == 1)))
text(c.r_time_elapsed(find(c.r_stage == 1)), h,"Ascent")
xline(c.r_time_elapsed(find(c.r_stage == 2)))
text(c.r_time_elapsed(find(c.r_stage == 2)), h,"Cruise")
xline(c.r_time_elapsed(find(c.r_stage == 3)))
text(c.r_time_elapsed(find(c.r_stage == 3)), h,"Descent")
xline(c.r_time_elapsed(find(c.r_stage == 4)))
text(c.r_time_elapsed(find(c.r_stage == 4)), h,"Landing")
% annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "pi value is " + pi)
%%
figure;
plot(c.r_time_elapsed,c.r_altitude,'.');
ylabel('Altitude [m]');
xlabel('Time [h]');
xline(c.r_time_elapsed(find(c.r_stage == 1)))
text(c.r_time_elapsed(find(c.r_stage == 1)), 0,"Ascent")
xline(c.r_time_elapsed(find(c.r_stage == 2)))
text(c.r_time_elapsed(find(c.r_stage == 2)), 0,"Cruise")
xline(c.r_time_elapsed(find(c.r_stage == 3)))
text(c.r_time_elapsed(find(c.r_stage == 3)), 0,"Descent")
xline(c.r_time_elapsed(find(c.r_stage == 4)))
text(c.r_time_elapsed(find(c.r_stage == 4)), 0,"Landing")
% annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "pi value is " + pi)

%%                                                                                                                                 
figure;
plot(c.r_time_elapsed,c.r_T);
ylabel('T');
xlabel('Time [h]');

figure;
plot(c.r_time_elapsed,c.r_P);
ylabel('P');
xlabel('Time [h]');

figure;
plot(c.r_time_elapsed,c.r_V);
ylabel('V');
xlabel('Time [h]');

figure;
plot(c.r_time_elapsed,c.r_eta);
ylabel('eta');
xlabel('Time [h]');

%%
figure;
plot(c.r_time_elapsed,c.r_kWh_avail);
ylabel('Energy Avaliable [kWh]');
xlabel('Time [h]');

%%
figure;
plot(c.r_time_elapsed,c.r_CL);
ylabel('CL');
xlabel('Time [h]');

figure;
plot(c.r_time_elapsed,c.r_CDi);
ylabel('CDi');
xlabel('Time [h]');

figure;
plot(c.r_time_elapsed,c.r_CD);
ylabel('CD');
xlabel('Time [h]');


%% compare to other aircraft - run mission

% set timestep, run mission
embraerLineage.tstep = tstep;
embraerLineage.mission_(m_length);
bae146.tstep = tstep;
bae146.mission_(m_length);
b747.tstep = tstep;
b747.mission_(m_length);
cub.tstep = tstep;
cub.mission_(m_length);
c130.tstep = tstep;
c130.mission_(m_length);

%% compare to other aircraft - plot results

% save open figures  later 
dest = '/home/ama/Documents/hybrid_plane_experts/CDR_figures/competition_comparison';
%dest = 'D:\MATLAB\AEM 4331\plots';
k = 0;

% compare properties
aircraft = [embraerLineage bae146 b747 cub c130 c];
a_names = [aircraft.Name];%['Embraer Lineage 1000' 'BAE 146' 'Boeing 747-8' 'J3C-65 Cub' 'C-130H']
N = categorical(a_names);
N = reordercats(N,a_names);

FigH = figure;
b = bar(N,[aircraft.time_elapsed]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Total Flight Time [h]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.flight_cost]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Total Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.c]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Chord Length [m]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.b]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Wing Span [m]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.M]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Total Mass [kg]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.Vcruise]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Cruise Speed [km/h]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
b = bar(N,[aircraft.kWh_avail]);
b.FaceColor = 'flat';
b.CData(end,:) = [.5 .8 .5];
ylabel('Remaining Energy after Flight [kWh')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

%% vary props of our aircraft (see prev. code)

v_aircraft = [];
n_var = 5;
k = 0;
dest = '/home/ama/Documents/hybrid_plane_experts/CDR_figures/properties_comparison';


% c
c_vary = c.c + linspace(-2,2,n_var)
for i = 1:n_var
    t = Aircraft;
    t.tstep = tstep;
    t.c = c_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.c], [v_aircraft.kWh_avail]);
hold on;
scatter(c.c,c.kWh_avail)
xlabel('c [m]')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.c], [v_aircraft.flight_cost]);
hold on;
scatter(c.c,c.flight_cost)
xlabel('c [m]')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

%% b
v_aircraft = [];
b_vary = c.b + linspace(-15,15,n_var);
for i = 1:n_var
    t = Aircraft;
    t.tstep = tstep;
    t.b = b_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.b], [v_aircraft.kWh_avail]);
hold on;
scatter(c.b,c.kWh_avail)
xlabel('b [m]')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.b], [v_aircraft.flight_cost]);
hold on;
scatter(c.b,c.flight_cost)
xlabel('b [m]')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

%% fratio
v_aircraft = [];
fratio_vary = c.fratio + linspace(-0.10,0.8,n_var);
for i = 1:n_var
    t = Aircraft;
    t.tstep = tstep;
    t.fratio = fratio_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.fratio], [v_aircraft.kWh_avail]);
hold on;
scatter(c.fratio,c.kWh_avail)
xlabel('Fuel Ratio')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.fratio], [v_aircraft.flight_cost]);
hold on;
scatter(c.fratio,c.flight_cost)
xlabel('Fuel Ratio')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;


%% dens_batt
v_aircraft = [];
dens_batt_vary = c.dens_batt + linspace(-0.5,1,n_var);
for i = 1:n_var
    t = Aircraft;
    t.tstep = tstep;
    t.dens_batt = dens_batt_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.dens_batt], [v_aircraft.kWh_avail]);
hold on;
scatter(c.dens_batt,c.kWh_avail)
xlabel('Battery Energy Density [kWh/kg]')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.dens_batt], [v_aircraft.flight_cost]);
hold on;
scatter(c.dens_batt,c.flight_cost)
xlabel('Battery Energy Density [kWh/kg]')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;


%% tstep
v_aircraft = [];
mult = 5;
tstep_vary = c.tstep*linspace(0.1,10,n_var*mult);
for i = 1:n_var*mult
    t = Aircraft;
    t.tstep = tstep_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.tstep], [v_aircraft.kWh_avail]);
hold on;
scatter(c.tstep,c.kWh_avail)
xlabel('Time step [h]')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.tstep], [v_aircraft.flight_cost]);
hold on;
scatter(c.tstep,c.flight_cost)
xlabel('Time step [h]')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;


%% 


%% d (prop diameter)
v_aircraft = [];
d_vary = c.d + linspace(-1.5,1.5,n_var);
for i = 1:n_var
    t = Aircraft;
    t.tstep = tstep;
    t.d = d_vary(i);
    t.mission_(m_length);
    v_aircraft = [v_aircraft t];
end

%%
FigH = figure;
plot([v_aircraft.d], [v_aircraft.kWh_avail]);
hold on;
scatter(c.d,c.kWh_avail)
xlabel('d (Prop. Diameter) [m]')
ylabel('Remaining Energy after Flight [kWh]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;

FigH = figure;
plot([v_aircraft.d], [v_aircraft.flight_cost]);
hold on;
scatter(c.d,c.flight_cost)
xlabel('d (Prop. Diameter) [m]')
ylabel('Flight Cost [$]')
FilePre = sprintf('%s-%d', 'MyFig', k);
savefig(FigH, fullfile(dest, [FilePre '.fig']))
print(FigH, fullfile(dest, [FilePre '.png']), '-r300', '-dpng')
k = k + 1;


%% save all open figures
% https://www.mathworks.com/matlabcentral/answers/182574-save-all-the-plots
% THIS IS BROKEN DON"T BOTHER THO
% FolderName = 'CDR_figures';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(FolderName, [FigName '.fig']));
% end
