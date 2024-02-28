clc


% 1 operating point
% 2 barometric pressure (mBar)
% 3 inlet pressure (kPa)
% 4 air flow pressure drop (Pa)
% 5 oil pressure (bar)
% 6 fuel pressure (bar)
% 7 throttle (%)
% 8 torque (Nm)
% 9 engine speed (RPM)
% 10 fuel flow rate (cm^3/min)
% 11 air temp (deg C)
% 12 fuel temp (deg C)
% 13 oil out temp (deg C)
% 14 exhaust temp (deg C)
% 15 plant in temp (deg C)
% 16 cool in temp (deg C)
% 17 cool out temp (deg C)
% 18 dyno out temp (deg C)
% 19 oil in temp (deg C)
% 20 cal water in temp (deg C)
% 21 cal water out temp (deg C)
% 22 cal exhaust in temp (deg C)
% 23 cal exhaust out temp (deg C)

Pp = petrol;
Pr = table2array(Pp);

Dp = diesel;
Dr = table2array(Dp);

shape = [6,23];
rhoa = 1.2;

%
% Petrol
vd_p = 1.2 * 10^(-3);
rv_p = 10.3;
rhof_p = 770;
lcvf_p = 44000;
afrst_p = 14.6;
clfm_p = 50.2 * 10^(-6);

P = zeros(5,9);

for t = 1:5
    P(t,1) = clfm_p * rhoa * Pr(t,4); % 1 air mass flowrate
    P(t,2) = rhof_p * Pr(t,10) * 10^(-6)/60; % 2 fuel mass flowrate
    P(t,3) = 0.001 * Pr(t,8) * 2 * pi * Pr(t,9) / 60; % 3 brake power (kW)
    P(t,4) = 0.001 * 4 * pi * Pr(t,8) / vd_p; % 4 brake mean effective pressure (kPa)
    P(t,5) = P(t,3) / (P(t,2) * lcvf_p); % 5 brake thermal efficiency
    P(t,6) = P(t,2) / P(t,3); % 6 brake specific fuel consumption
    P(t,7) = 120 * clfm_p * Pr(t,4) / (vd_p * Pr(t,9)); % 7 volumetric efficiency
    P(t,8) = P(t,1) / P(t,2); % 8 air/fuel ratio
    P(t,9) = 1 - rv_p^(-0.4); % 9 theoretical efficiency
end
%

%
% Diesel
vd_d = 1.9 * 10^(-3);
rv_d = 22;
rhof_d = 840;
lcvf_d = 42500;
afrst_d = 14.5;
clfm_d = 55.9 * 10^(-6);

D = zeros(5,9);

for t = 1:5
    D(t,1) = clfm_p * rhoa * Dr(t,4); % 1 air mass flowrate
    D(t,2) = rhof_p * Dr(t,10) * 10^(-6)/60; % 2 fuel mass flowrate
    D(t,3) = 0.001 * Dr(t,8) * 2 * pi * Dr(t,9) / 60; % 3 brake power (kW)
    D(t,4) = 0.001 * 4 * pi * Dr(t,8) / vd_p; % 4 brake mean effective pressure (kPa)
    D(t,5) = D(t,3) / (D(t,2) * lcvf_p); % 5 brake thermal efficiency
    D(t,6) = D(t,2) / D(t,3); % 6 brake specific fuel consumption
    D(t,7) = 120 * clfm_p * Dr(t,4) / (vd_p * Dr(t,9)); % 7 volumetric efficiency
    D(t,8) = D(t,1) / D(t,2); % 8 air/fuel ratio
    D(t,9) = (Dr(t,14)/Dr(t,11))^(1/1.4); % 9 cutoff ratio
    D(t,10) = 1 - rv_p^(-0.4); % 10 theoretical efficiency
end
%



% Plotting
figure
hold on
box on
grid on
xlabel('Engine speed (RPM)');

plot(Pr(:,9),P(:,4), '--');
scatter(Pr(:,9),P(:,4), '*');
ylabel('Bmep (kPa)');
