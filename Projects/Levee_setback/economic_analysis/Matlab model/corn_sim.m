%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program   : corn_sim.m
% Programmer: Yusuke Kuwayama
% Date      : June 7, 2023
% Ref       : Simulate profits from corn production given a groundwater and
%             surface water irrigation schedule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load data %%%
load("K_C.mat") % Load daily crop coefficients
load("ET_O.mat") % Load daily reference ET (in)
load("PRECIP.mat") % Load daily precipitation (in)
load("RO.mat") % Load daily runoff (in)
load("CR.mat") % Load daily capillary rise (in)
load("DP.mat") % Load daily deep percolation (in)
load("I_GW.mat") % Load daily groundwater irrigation (in)
load("I_SW.mat") % Load daily surface water irrigation (in)
load("K_Y.mat") % Load daily yield response factors
load("X.mat") % Load daily groundwater depths (ft)

%%% Set parameter values %%%
season_start = datetime(2022,4,30); % Start date for growing season
season_end = datetime(2022,9,8); % End date for growing season
p_table22 = 0.55; % Soil water depletion fraction for no stress
theta_fc = 0.29; % Water content at field capacity
theta_wp = 0.15; % Water content at wilting point
zr = 53.15; % Rooting depth (in)
d_ini = 0; % Initial root zone depletion (in)
y_max = 6; % Maximum expected crop yield (tons/acre)
phi = 0.13; % Energy requirement to raise a unit of water by a unit of vertical distance (kWh/acre-in/ft)
p_c = 210; % Crop price ($/ton)
p_sw = 2; % Surface water charges and fees ($/acre-in)
p_e = 0.17; % Cost of energy for groundwater pumping ($/kWh)
p_o = 956; % Variable operating costs per acre, excluding irrigation costs ($/acre)
a = 1; % Parcel area

%%% Create vectors of growing season dates %%%
DATE = transpose(season_start:season_end); % Create vector of dates for growing season
SEASON_DAY = transpose(0:length(DATE)-1); % Create vector of days elapsed in growing season

%%% Calculate daily root zone depletion (in) %%%
taw = (theta_fc - theta_wp)*zr; % Calculate total available water in the root zone
ET_C = K_C.*ET_O; % Calculate daily crop ET under normal conditions
P = p_table22 + 0.04*((5-(25.4*ET_C))); % Calculate adjusted daily soil water depletion fraction for no stress
RAW = taw*P; % Calculate readily available water in the root zone
D = zeros(length(DATE),1); % Create empty vector for daily root zone depletion
D(1,1) = min(d_ini + ET_C(1,1) - PRECIP(1,1) + RO(1,1) + CR(1,1) + DP(1,1) - I_GW(1,1) - I_SW(1,1),taw); % Calculate root zone depletion for first day of the season (cannot exceed total available water)
for i = 2:length(DATE)
    D(i,1) = min(D(i-1,1) + ET_C(i,1) - PRECIP(i,1) + RO(i,1) + CR(i,1) + DP(i,1) - I_GW(i,1) - I_SW(i,1),taw);
end % Calculate root zone depletion for remaining days of the season (cannot exceed total available water)

%%% Calculate daily water stress coefficient %%%
K_S = zeros(length(DATE),1); % Create empty vector for daily water stress coefficient
for i = 1:length(DATE)
    if D(i,1) > RAW(i,1)
        K_S(i,1) = (taw - D(i,1))/((1 - P(i,1))*taw);
    else
        K_S(i,1) = 1;
    end
end % Calculate daily water stress coefficient (equals 1 if there is no soil water stress)

%%% Calculate daily crop outcomes %%%
ET_CADJ = K_S.*K_C.*ET_O; % Calculate daily crop ET with soil water stress
Y_A = y_max*(ones(length(DATE),1) - K_Y.*(ones(length(DATE),1) - (ET_CADJ./ET_C))); % Calculate actual yield

%%% Calculate economic outcomes %%%
C_GW = p_e*phi*(X.*I_GW); % Calculate daily groundwater pumping costs ($/acre)
y_atot = mean(Y_A); % Calculate total yield for the season (tons/acre)
i_gwtot = sum(I_GW); % Calculate total groundwater irrigation for the season (in)
i_swtot = sum(I_SW); % Calculate total surface water irrigation for the season (in)
c_gwtot = sum(C_GW); % Calcualte total groundwater pumping costs for the season ($/acre)
c_swtot = p_sw*i_swtot; % Calcualte total surface water costs for the season ($/acre)
pi = a*(p_c*y_atot - c_gwtot - c_swtot - p_o); % Calculate profit  ($)

subplot(3,1,1)
plot(DATE,ET_O)
title('Daily reference ET (in)')
ylim([0 0.5])
subplot(3,1,2)
plot(DATE,ET_C)
title('Daily crop ET under normal conditions')
ylim([0 0.5])
subplot(3,1,3)
plot(DATE,ET_CADJ)
title('Daily crop ET with soil water stress')
ylim([0 0.5])