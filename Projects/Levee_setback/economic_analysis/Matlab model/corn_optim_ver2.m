%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program   : corn_optim_ver2.m
% Programmer: Yusuke Kuwayama
% Date      : June 9, 2023
% Ref       : Find unconstrained and constrained optimal irrigation
%             schedule for corn production for surface water and
%             groundwater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set global variables %%%
global DATE n_irr IRR_DAYS d_ini ET_C taw RAW P K_C ET_O y_max p_e phi p_sw a p_c p_o

%%% Load data %%%
load("K_C.mat") % Load daily crop coefficients
load("ET_O.mat") % Load daily reference ET (in)

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
gap_irr = 14; % Number of days between irrigations

% adjust here for sensitivity testing
I_WTOTCON = [5 10]'; % Total surface water and groundwater available during the season (in)
%%% Create vectors of growing season dates %%%
DATE = transpose(season_start:season_end); % Create vector of dates for growing season
n_irr = floor(length(DATE)/gap_irr) + 1; % Calculate number of irrigations
IRR_DAYS = (0:gap_irr:(n_irr*gap_irr-1))'; % Calculate days on which irrigation takes place

%%% Initial calculations %%%
taw = (theta_fc - theta_wp)*zr; % Calculate total available water in the root zone
ET_C = K_C.*ET_O; % Calculate daily crop ET under normal conditions
P = p_table22 + 0.04*((5-(25.4*ET_C))); % Calculate adjusted daily soil water depletion fraction for no stress
RAW = taw*P; % Calculate readily available water in the root zone

%%% Calculate unconstrained optimal irrigation schedule %%%
I_WMAX0UNCON = zeros(1,2*n_irr); % Initial irrigation values for optimization
AUNCON = []; % No inequality constraints
I_WTOTUNCON = []; % No inequality constraints
AEQUNCON = []; % No equality contraints
BEQUNCON = []; % No equality contraints
I_WMAXLBUNCON = zeros(1,2*n_irr); % Irrigation cannot be negative
[I_WMAXUNCON,PI_MAXUNCON] = fmincon(@corn_profit_ver2, I_WMAX0UNCON, AUNCON, I_WTOTUNCON, AEQUNCON, BEQUNCON, I_WMAXLBUNCON); % Calculate optimal irrigations using profit function
PI_MAXUNCON = (-1)*PI_MAXUNCON; % Make profit positive

%%% Calculate optimal irrigation schedule %%%
I_WMAX0CON = zeros(1,2*n_irr); % Initial irrigation values for optimization
% Coefficients for inequality constraints (first n_irr columns are for surface water; second n_irr columns are for groundwater)
ACON = zeros(1,n_irr);
ACON(1,1:n_irr) = ones(1,n_irr); 
ACON(2,(n_irr+1):(2*n_irr)) = ones(1,n_irr);
AEQCON = []; % No equality contraints
BEQCON = []; % No equality contraints
I_WMAXLBCON = zeros(1,2*n_irr); % Irrigation cannot be negative
[I_WMAXCON,PI_MAXCON] = fmincon(@corn_profit_ver2, I_WMAX0CON, ACON, I_WTOTCON, AEQCON, BEQCON, I_WMAXLBCON); % Calculate optimal irrigations using profit function
PI_MAXCON = (-1)*PI_MAXCON; % Make profit positive

%% Plot irrigation during the season
IRR_DATES = NaT(n_irr,1);
for i = 1:n_irr
    IRR_DATES(i,1) = DATE(IRR_DAYS(i,1)+1,1);
end
subplot(2,1,1)
bar(IRR_DATES,[I_WMAXUNCON(1:n_irr) ; I_WMAXUNCON((n_irr+1):(2*n_irr))]')
title("Unconstrained optimal water use")
ylim([0 ceil(max([I_WMAXUNCON(1:n_irr) I_WMAXUNCON((n_irr+1):(2*n_irr))]))])
ylabel('Irrigation (in)')
legend('surface water','groundwater')
subplot(2,1,2)
bar(IRR_DATES,[I_WMAXCON(1:n_irr) ; I_WMAXCON((n_irr+1):(2*n_irr))]')
title("Constrained optimal water use (w^{SW}=" + I_WTOTCON(1,1) + ", w^{GW}=" + I_WTOTCON(2,1) + ")")
ylim([0 ceil(max([I_WMAXUNCON(1:n_irr) I_WMAXUNCON((n_irr+1):(2*n_irr))]))])
ylabel('Irrigation (in)')
legend('surface water','groundwater')