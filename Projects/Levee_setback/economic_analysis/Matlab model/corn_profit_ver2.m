%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program   : corn_profit_ver2.m
% Programmer: Yusuke Kuwayama
% Date      : June 9, 2023
% Ref       : Profit function for corn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pi = corn_profit_ver2(I_WLEVELS)

%%% Set global variables %%%
global DATE n_irr IRR_DAYS d_ini ET_C taw RAW P K_C ET_O y_max p_e phi p_sw a p_c p_o

%%% Load data %%%
load("PRECIP.mat"); % Load daily precipitation (in)
load("RO.mat") % Load daily runoff (in)
load("CR.mat") % Load daily capillary rise (in)
load("DP.mat") % Load daily deep percolation (in)
load("K_Y.mat") % Load daily yield response factors
load("X.mat") % Load daily groundwater depths (ft)

%%% Define which days in the season are irrigation days %%%
I_SW = zeros(length(DATE),1);
I_GW = zeros(length(DATE),1);
for i = 1:n_irr
    I_SW(IRR_DAYS(i,1)+1,1) = I_WLEVELS(i);
    I_GW(IRR_DAYS(i,1)+1,1) = I_WLEVELS(i+n_irr);
end

%%% Calculate daily root zone depletion (in) %%%
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
y_atot = mean(Y_A); % Calculate total yield for the season (tons/acre)
c_gwtot = sum(p_e*phi*(X.*I_GW)); % Calcualte total groundwater pumping costs for the season ($/acre)
c_swtot = p_sw*sum(I_SW); % Calcualte total surface water costs for the season ($/acre)
pi = -(a*(p_c*y_atot - c_gwtot - c_swtot - p_o)); % Calculate profit ($)