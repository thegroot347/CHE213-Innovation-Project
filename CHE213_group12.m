%% Finding the calibration parameters
clc; close all; clear;
% Data
Feed_RI = [1.36120 1.36399 1.36408 1.36958];
hexane_v_v_percentage = [0.2 0.4 0.6 0.8];
% Linear regression (y = mx + c)
p = polyfit(Feed_RI, hexane_v_v_percentage, 1);
slope = p(1);
intercept = p(2);
% Evaluating fitted line for plotting
y_fit = polyval(p, Feed_RI);
% Calculating R-squared (R^2)
y_resid = hexane_v_v_percentage - y_fit;
SSresid = sum(y_resid.^2);
SStotal = (length(hexane_v_v_percentage)-1) * var(hexane_v_v_percentage);
R2 = 1 - SSresid/SStotal;

% Plotting Calibration
figure('Color', 'w');
movegui(gcf, 'center'); % Centers the figure on screen
hold on;
plot(Feed_RI, hexane_v_v_percentage, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5); 
plot(Feed_RI, y_fit, '-r', 'LineWidth', 2);                             
xlabel('Refractive Index');
ylabel('Volume Fraction (v/v)');
title(['Calibration Curve (R^2 = ', num2str(R2, '%.4f'), ')']);
legend('Experimental Data', ['Fit: y = ', num2str(slope, '%.2f'), 'x ', num2str(intercept, '%.2f')], 'Location', 'NorthWest');
grid on;
hold off;

%% Properties and Functions
rho_hex = 659.0;    M_hex = 86.18;      % Density at 25°C (kg/m^3), MW (g/mol)
rho_eth = 785.1;    M_eth = 46.07;      % Density at 25°C (kg/m^3), MW (g/mol)
k_hex = rho_hex / M_hex;
k_eth = rho_eth / M_eth;
% Function expects v_frac as a decimal (0 to 1)
mol_fraction_hexane = @(v_frac) (v_frac * k_hex) / (v_frac * k_hex + (1 - v_frac) * k_eth);

%% Sample 1 (Hexane: 20% v/v nominal)
RI_top_1 = [1.36949 1.36983 1.36973 1.36048];
RI_bot_1 = [1.36093 1.36088 1.36048 1.36090];
RI_top_1_avg = mean(RI_top_1);
RI_bot_1_avg = mean(RI_bot_1);
v_per_bot1 = slope * RI_bot_1_avg + intercept;
v_per_top1 = slope * RI_top_1_avg + intercept;
x_hexane_1 = mol_fraction_hexane(v_per_bot1); 
y_hexane_1 = mol_fraction_hexane(v_per_top1); 

%% Sample 2 (Hexane: 40% v/v nominal)
RI_top_2 = [1.37151 1.37139 1.37142 1.37146];
RI_bot_2 = [1.36382 1.36365 1.36330 1.36341];
RI_top_2_avg = mean(RI_top_2);
RI_bot_2_avg = mean(RI_bot_2);
v_per_bot2 = slope * RI_bot_2_avg + intercept;
v_per_top2 = slope * RI_top_2_avg + intercept;
x_hexane_2 = mol_fraction_hexane(v_per_bot2); 
y_hexane_2 = mol_fraction_hexane(v_per_top2); 

%% Sample 3 (Hexane: 60% v/v nominal)
RI_top_3 = [1.37087 1.37089 1.37080 1.37096];
RI_bot_3 = [1.36730 1.36714 1.36746 1.36739];
RI_top_3_avg = mean(RI_top_3);
RI_bot_3_avg = mean(RI_bot_3);
v_per_bot3 = slope * RI_bot_3_avg + intercept;
v_per_top3 = slope * RI_top_3_avg + intercept;
x_hexane_3 = mol_fraction_hexane(v_per_bot3); 
y_hexane_3 = mol_fraction_hexane(v_per_top3); 

%% Sample 4 (Hexane: 80% v/v nominal)
RI_top_4 = [1.37157 1.37162 1.37157 1.37171];
RI_bot_4 = [1.37148 1.37157 1.37148 1.37152];
RI_top_4_avg = mean(RI_top_4);
RI_bot_4_avg = mean(RI_bot_4);
v_per_bot4 = slope * RI_bot_4_avg + intercept;
v_per_top4 = slope * RI_top_4_avg + intercept;
x_hexane_4 = mol_fraction_hexane(v_per_bot4); 
y_hexane_4 = mol_fraction_hexane(v_per_top4); 

%% VLE Modeling: Experiment vs Wilson, NRTL, and UNIQUAC
%% 1. Pure Component & Model Parameters
R = 1.987; % Gas constant in cal/(mol*K)
P_total = 760; % System pressure in mmHg
% Antoine Constants
A1 = 6.87776; B1 = 1171.53; C1 = 224.366;
A2 = 8.04494; B2 = 1554.3;  C2 = 222.65;

% Wilson Parameters
V1 = 131.61; V2 = 58.68;           
W_a12 = 437.98; W_a21 = 1438.0;    

% NRTL Parameters
alpha = 0.47;                        
N_dg12 = 1036.0; N_dg21 = 1095.0;  

% UNIQUAC Parameters
r1 = 4.4998; q1 = 3.856;   
r2 = 2.1055; q2 = 1.972;   
z = 10; 
U_a12 = 328.34; U_a21 = 339.74; 

%% 2. Calculating Theoretical Curves (x-y Generation)
x_fit = linspace(0.001, 0.999, 100);
y_wilson = zeros(size(x_fit));
y_nrtl   = zeros(size(x_fit));
y_uniquac= zeros(size(x_fit));

for i = 1:length(x_fit)
    x1 = x_fit(i); x2 = 1 - x1;
    
    % WILSON
    T = x1*68.7 + x2*78.3; 
    for iter = 1:5 
        T_K = T + 273.15;
        L12 = (V2/V1) * exp(-W_a12 / (R * T_K));
        L21 = (V1/V2) * exp(-W_a21 / (R * T_K));
        gamma1 = exp(-log(x1 + x2*L12) + x2 * (L12/(x1 + x2*L12) - L21/(x2 + x1*L21)));
        Psat1 = 10^(A1 - B1/(T + C1));
        y_wilson(i) = (x1 * gamma1 * Psat1) / P_total;
        T = T + 0.1; 
    end
    
    % NRTL
    T = x1*68.7 + x2*78.3;
    for iter = 1:5
        T_K = T + 273.15;
        tau12 = N_dg12 / (R * T_K); tau21 = N_dg21 / (R * T_K);
        G12 = exp(-alpha * tau12); G21 = exp(-alpha * tau21);
        gamma1 = exp(x2^2 * (tau21 * (G21 / (x1 + x2*G21))^2 + (tau12 * G12) / (x2 + x1*G12)^2));
        Psat1 = 10^(A1 - B1/(T + C1));
        y_nrtl(i) = (x1 * gamma1 * Psat1) / P_total;
        T = T + 0.1;
    end
    
    % UNIQUAC
    T = x1*68.7 + x2*78.3;
    for iter = 1:5
        T_K = T + 273.15;
        phi1 = (x1*r1) / (x1*r1 + x2*r2); phi2 = (x2*r2) / (x1*r1 + x2*r2);
        th1  = (x1*q1) / (x1*q1 + x2*q2); th2  = (x2*q2) / (x1*q1 + x2*q2);
        l1 = (z/2)*(r1 - q1) - (r1 - 1); l2 = (z/2)*(r2 - q2) - (r2 - 1);
        ln_gC1 = log(phi1/x1) + (z/2)*q1*log(th1/phi1) + l1 - (phi1/x1)*(x1*l1 + x2*l2);
        tau12 = exp(-U_a12 / (R * T_K)); tau21 = exp(-U_a21 / (R * T_K));
        ln_gR1 = q1 * (1 - log(th1 + th2*tau21) - (th1 / (th1 + th2*tau21)) - (th2*tau12 / (th2 + th1*tau12)));
        gamma1 = exp(ln_gC1 + ln_gR1);
        Psat1 = 10^(A1 - B1/(T + C1));
        y_uniquac(i) = (x1 * gamma1 * Psat1) / P_total;
        T = T + 0.1;
    end
end

%% 3. Experimental Data & Constrained Polynomial Fit
x_exp = [x_hexane_1 x_hexane_2 x_hexane_3 x_hexane_4]; 
y_exp = [y_hexane_1 y_hexane_2 y_hexane_3 y_hexane_4];

% Experimental Temperatures 
T_liq_exp = [67.9, 59.7, 59.1, 56.8]; % T_Bottom
T_vap_exp = [55.7, 55.5, 54.9, 54.8]; % T_Top

% Calculating Experimental Activity Coefficients
gamma1_exp = zeros(1,4);
gamma2_exp = zeros(1,4);
for i = 1:4
    T = T_liq_exp(i);
    Psat1 = 10^(A1 - B1/(T + C1));
    Psat2 = 10^(A2 - B2/(T + C2));
    gamma1_exp(i) = (y_exp(i) * P_total) / (x_exp(i) * Psat1);
    gamma2_exp(i) = ((1 - y_exp(i)) * P_total) / ((1 - x_exp(i)) * Psat2);
end

poly_coeffs = [-2.7952, 10.2787, -12.1817, 5.6982, 0.0000];
y_poly_fit = polyval(poly_coeffs, linspace(0, 1, 100));

%% 4. Final Combined Plot (x-y Diagram)
figure('Color', 'w', 'Position', [0, 0, 800, 600]);
movegui(gcf, 'center'); % Centers the figure on screen
hold on;
plot(x_fit, y_wilson,  'b-', 'LineWidth', 2, 'DisplayName', 'Wilson Model');
plot(x_fit, y_nrtl,    'r--', 'LineWidth', 2, 'DisplayName', 'NRTL Model');
plot(x_fit, y_uniquac, 'k-.', 'LineWidth', 2, 'DisplayName', 'UNIQUAC Model');
plot(linspace(0, 1, 100), y_poly_fit, 'm:', 'LineWidth', 2.5, 'DisplayName', '4th-Order Poly Fit');
plot(x_exp, y_exp, '*', 'Color', [0 0.5 0], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Experimental Data');
plot([0 1], [0 1], 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility', 'off'); 

xlabel('x_{Hexane} (Liquid Mole Fraction)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y_{Hexane} (Vapor Mole Fraction)', 'FontSize', 12, 'FontWeight', 'bold');
title('Hexane-Ethanol VLE (x-y Diagram)', 'FontSize', 14);
grid on; legend('Location', 'southeast', 'FontSize', 11);
xlim([0 1]); ylim([0 1]); hold off;

%% 5. Exact T-xy Diagram and Activity Coefficient Calculations 
T_bub_wilson = zeros(size(x_fit)); y_dew_wilson = zeros(size(x_fit));
T_bub_nrtl   = zeros(size(x_fit)); y_dew_nrtl   = zeros(size(x_fit));
T_bub_uniquac= zeros(size(x_fit)); y_dew_uniquac= zeros(size(x_fit));

gamma1_wilson = zeros(size(x_fit)); gamma2_wilson = zeros(size(x_fit));
gamma1_nrtl   = zeros(size(x_fit)); gamma2_nrtl   = zeros(size(x_fit));
gamma1_uniquac= zeros(size(x_fit)); gamma2_uniquac= zeros(size(x_fit));

for i = 1:length(x_fit)
    x1 = x_fit(i); x2 = 1 - x1;
    
    % Wilson
    T = x1*68.7 + x2*78.3; 
    for iter=1:20
        T_K = T + 273.15;
        L12 = (V2/V1) * exp(-W_a12 / (R * T_K));
        L21 = (V1/V2) * exp(-W_a21 / (R * T_K));
        g1 = exp(-log(x1 + x2*L12) + x2 * (L12/(x1 + x2*L12) - L21/(x2 + x1*L21)));
        g2 = exp(-log(x2 + x1*L21) + x1 * (L21/(x2 + x1*L21) - L12/(x1 + x2*L12)));
        P1 = 10^(A1 - B1/(T + C1)); P2 = 10^(A2 - B2/(T + C2));
        err = (x1*g1*P1 + x2*g2*P2) - P_total;
        if abs(err) < 0.01, break; end
        dP_dT = x1*g1*(P1*B1*log(10)/(T+C1)^2) + x2*g2*(P2*B2*log(10)/(T+C2)^2);
        T = T - err / dP_dT; 
    end
    T_bub_wilson(i) = T; y_dew_wilson(i) = (x1 * g1 * P1) / P_total;
    gamma1_wilson(i) = g1; gamma2_wilson(i) = g2;

    % NRTL
    T = x1*68.7 + x2*78.3; 
    for iter=1:20
        T_K = T + 273.15;
        tau12 = N_dg12 / (R * T_K); tau21 = N_dg21 / (R * T_K);
        G12 = exp(-alpha * tau12); G21 = exp(-alpha * tau21);
        g1 = exp(x2^2 * (tau21 * (G21 / (x1 + x2*G21))^2 + (tau12 * G12) / (x2 + x1*G12)^2));
        g2 = exp(x1^2 * (tau12 * (G12 / (x2 + x1*G12))^2 + (tau21 * G21) / (x1 + x2*G21)^2));
        P1 = 10^(A1 - B1/(T + C1)); P2 = 10^(A2 - B2/(T + C2));
        err = (x1*g1*P1 + x2*g2*P2) - P_total;
        if abs(err) < 0.01, break; end
        dP_dT = x1*g1*(P1*B1*log(10)/(T+C1)^2) + x2*g2*(P2*B2*log(10)/(T+C2)^2);
        T = T - err / dP_dT; 
    end
    T_bub_nrtl(i) = T; y_dew_nrtl(i) = (x1 * g1 * P1) / P_total;
    gamma1_nrtl(i) = g1; gamma2_nrtl(i) = g2;
    
    % UNIQUAC
    T = x1*68.7 + x2*78.3; 
    for iter=1:20
        T_K = T + 273.15;
        phi1 = (x1*r1)/(x1*r1 + x2*r2); phi2 = (x2*r2)/(x1*r1 + x2*r2);
        th1  = (x1*q1)/(x1*q1 + x2*q2); th2  = (x2*q2)/(x1*q1 + x2*q2);
        l1 = (z/2)*(r1 - q1) - (r1 - 1); l2 = (z/2)*(r2 - q2) - (r2 - 1);
        ln_gC1 = log(phi1/x1) + (z/2)*q1*log(th1/phi1) + l1 - (phi1/x1)*(x1*l1 + x2*l2);
        ln_gC2 = log(phi2/x2) + (z/2)*q2*log(th2/phi2) + l2 - (phi2/x2)*(x1*l1 + x2*l2);
        tau12 = exp(-U_a12 / (R * T_K)); tau21 = exp(-U_a21 / (R * T_K));
        ln_gR1 = q1 * (1 - log(th1 + th2*tau21) - (th1/(th1 + th2*tau21)) - (th2*tau12/(th2 + th1*tau12)));
        ln_gR2 = q2 * (1 - log(th2 + th1*tau12) - (th2/(th2 + th1*tau12)) - (th1*tau21/(th1 + th2*tau21)));
        g1 = exp(ln_gC1 + ln_gR1); g2 = exp(ln_gC2 + g2); 
        g2 = exp(ln_gC2 + ln_gR2);
        P1 = 10^(A1 - B1/(T + C1)); P2 = 10^(A2 - B2/(T + C2));
        err = (x1*g1*P1 + x2*g2*P2) - P_total;
        if abs(err) < 0.01, break; end
        dP_dT = x1*g1*(P1*B1*log(10)/(T+C1)^2) + x2*g2*(P2*B2*log(10)/(T+C2)^2);
        T = T - err / dP_dT; 
    end
    T_bub_uniquac(i) = T; y_dew_uniquac(i) = (x1 * g1 * P1) / P_total;
    gamma1_uniquac(i) = g1; gamma2_uniquac(i) = g2;
end

%% 6. Plotting T-xy Diagrams 
figure('Color', 'w', 'Position', [0, 0, 1400, 450]);
movegui(gcf, 'center'); % Centers the figure on screen

subplot(1, 3, 1); hold on;
plot(x_fit, T_bub_wilson, '-r', 'LineWidth', 2, 'DisplayName', 'Wilson Bubble Curve');
plot(y_dew_wilson, T_bub_wilson, '-b', 'LineWidth', 2, 'DisplayName', 'Wilson Dew Curve');
plot(x_exp, T_liq_exp, 'rs', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Exp. Bubble Pt (x)');
plot(y_exp, T_vap_exp, 'bd', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Exp. Dew Pt (y)');
xlabel('x, y (Hexane Mole Fraction)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 10, 'FontWeight', 'bold');
title('Wilson', 'FontSize', 12); grid on; legend('Location', 'best', 'FontSize', 9); 
xlim([0 1]); ylim([50 80]); hold off; 

subplot(1, 3, 2); hold on;
plot(x_fit, T_bub_nrtl, '-r', 'LineWidth', 2, 'HandleVisibility', 'off'); 
plot(y_dew_nrtl, T_bub_nrtl, '-b', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_exp, T_liq_exp, 'rs', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(y_exp, T_vap_exp, 'bd', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('x, y (Hexane Mole Fraction)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 10, 'FontWeight', 'bold');
title('NRTL', 'FontSize', 12); grid on; 
xlim([0 1]); ylim([50 80]); hold off; 

subplot(1, 3, 3); hold on;
plot(x_fit, T_bub_uniquac, '-r', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(y_dew_uniquac, T_bub_uniquac, '-b', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_exp, T_liq_exp, 'rs', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(y_exp, T_vap_exp, 'bd', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('x, y (Hexane Mole Fraction)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 10, 'FontWeight', 'bold');
title('UNIQUAC', 'FontSize', 12); grid on; 
xlim([0 1]); ylim([50 80]); hold off; 

sgtitle('Hexane-Ethanol T-xy Diagrams: Phase-Consistent Color Scheme', 'FontSize', 14, 'FontWeight', 'bold');

%% 7. Plotting Activity Coefficients vs Composition
figure('Color', 'w', 'Position', [0, 0, 1400, 450]);
movegui(gcf, 'center'); % Centers the figure on screen

subplot(1, 3, 1); hold on;
plot(x_fit, gamma1_wilson, '-b', 'LineWidth', 2, 'DisplayName', '\gamma_{Hex} (Model)');
plot(x_fit, gamma2_wilson, '-r', 'LineWidth', 2, 'DisplayName', '\gamma_{Eth} (Model)');
plot(x_exp, gamma1_exp, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'Exp \gamma_{Hex}');
plot(x_exp, gamma2_exp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Exp \gamma_{Eth}');
plot([0 1], [1 1], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
xlabel('x_{Hexane}', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Activity Coefficient (\gamma)', 'FontSize', 10, 'FontWeight', 'bold');
title('Wilson', 'FontSize', 12); grid on; legend('Location', 'best', 'FontSize', 9); 
xlim([0 1]); ylim([0 8]); hold off; 

subplot(1, 3, 2); hold on;
plot(x_fit, gamma1_nrtl, '-b', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_fit, gamma2_nrtl, '-r', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_exp, gamma1_exp, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(x_exp, gamma2_exp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot([0 1], [1 1], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
xlabel('x_{Hexane}', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Activity Coefficient (\gamma)', 'FontSize', 10, 'FontWeight', 'bold');
title('NRTL', 'FontSize', 12); grid on; 
xlim([0 1]); ylim([0 8]); hold off; 

subplot(1, 3, 3); hold on;
plot(x_fit, gamma1_uniquac, '-b', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_fit, gamma2_uniquac, '-r', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(x_exp, gamma1_exp, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(x_exp, gamma2_exp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot([0 1], [1 1], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
xlabel('x_{Hexane}', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Activity Coefficient (\gamma)', 'FontSize', 10, 'FontWeight', 'bold');
title('UNIQUAC', 'FontSize', 12); grid on; 
xlim([0 1]); ylim([0 8]); hold off; 

sgtitle('Activity Coefficients: Comparison of Local Composition Models', 'FontSize', 14, 'FontWeight', 'bold');