%% Thin-Film TMM Optical Simulation for spectral selective absorber/emitter

% Author: Nima Talebzadeh

% Description: This project is about Transfer Matrix Method (TMM) Optical Simulation for spectral selective absorber/emitter in thermo-photovoltaic (TPV) systems

% Calculates reflectance, transmittance, and emissivity of multilayer photonic structures.

% Transfer Matrix Method with thickness sweep - OPTIMIZED VERSION (FIXED)
clc;
clear;

% Example user inputs
tungstenFile = 'hfn.xlsx';
hfnFile = 'sic.xlsx';
thicknessTungsten = 10000; % Example thickness in nm

% Load data from Excel files
tungstenData = readtable(tungstenFile);
hfnData = readtable(hfnFile);

% Extract wavelength and complex refractive indices (assumed pre-aligned)
wavelength = tungstenData{:,1}; % Shared wavelength axis
n_tungsten = tungstenData{:,2};
k_tungsten = tungstenData{:,3};
n_hfn = hfnData{:,2};
k_hfn = hfnData{:,3};

% OPTIMIZATION 1: Reduce number of points in your sweeps
% Define thickness range for HfN layer with fewer points
fine_thickness_range = linspace(0.01, 100, 300); % Reduced from 20
coarse_thickness_range = linspace(101, 1000, 400); % Reduced from 400
thickness_range = [fine_thickness_range, coarse_thickness_range]; 

% Temperature range with fewer points
temperature_range = linspace(1000, 2500, 4); % Reduced from 200

% Initialize R, T, A arrays
reflection = zeros(length(thickness_range), length(wavelength));
transmission = zeros(length(thickness_range), length(wavelength));
absorption = zeros(length(thickness_range), length(wavelength));

% Create waitbar for progress tracking
waitbar_handle = waitbar(0, 'Starting calculations...', 'Name', 'Transfer Matrix Method Progress');

% Constants - renamed to avoid variable conflicts
planck_const = 6.626e-34; 
speed_light = 3e8; 
boltzmann_const = 1.3806e-23;

% Convert wavelength to meters once outside of loops
wavelength_m = wavelength * 1e-9;

% Pre-compute complex refractive indices to reduce repeated calculations
n_hfn_complex = n_hfn + 1i * k_hfn;
n_w_complex = n_tungsten + 1i * k_tungsten;

% Sweep thickness and calculate R, T, A
total_iterations = length(thickness_range);
for t_idx = 1:total_iterations
    thicknessHfN = thickness_range(t_idx);
    
    % Update the waitbar
    waitbar(t_idx/total_iterations, waitbar_handle, sprintf('Processing thickness: %.2f nm (%.1f%%)', ...
        thicknessHfN, 100*t_idx/total_iterations));
    
    % Pre-compute deltas for current thickness
    delta1 = 2*pi*n_hfn_complex .* thicknessHfN ./ wavelength;
    delta2 = 2*pi*n_w_complex .* thicknessTungsten ./ wavelength;
    
    for idx = 1:length(wavelength)
        q1 = n_hfn_complex(idx);
        q2 = n_w_complex(idx);


       M_hfn = [cos(delta1(idx)), -1i*sin(delta1(idx))/q1;
                 -1i*q1*sin(delta1(idx)), cos(delta1(idx))];

        M_w = [cos(delta2(idx)), -1i*sin(delta2(idx))/q2;
               -1i*q2*sin(delta2(idx)), cos(delta2(idx))];

        M = M_hfn * M_w;

        r = (M(1,1) + M(1,2) - M(2,1) - M(2,2)) / ...
            (M(1,1) + M(1,2) + M(2,1) + M(2,2));
        t = 2 / (M(1,1) + M(1,2) + M(2,1) + M(2,2));

        R = abs(r)^2;
        T = abs(t)^2;
        A = 1 - R - T;

        reflection(t_idx, idx) = R;
        transmission(t_idx, idx) = T;
        absorption(t_idx, idx) = A;
    end
end

% Close the waitbar
close(waitbar_handle);

% Initialize arrays
total_emissive_power_structure = zeros(length(thickness_range), length(temperature_range));
total_emissive_power_structure_BB = zeros(length(thickness_range), length(temperature_range));
Loss_Rate = zeros(length(thickness_range), length(temperature_range));

% Create waitbar for thermal calculations
thermal_waitbar = waitbar(0, 'Starting thermal calculations...', 'Name', 'Thermal Calculations Progress');

% Pre-compute wavelength differences for integration
dlambda = zeros(size(wavelength_m));
for w_idx = 1:length(wavelength_m)-1
    dlambda(w_idx) = wavelength_m(w_idx+1) - wavelength_m(w_idx);
end
dlambda(end) = dlambda(end-1); % Last point uses same interval as previous

total_calcs = length(thickness_range) * length(temperature_range);
calc_counter = 0;


% Extract thickness for each wavelength where absorption is maximum
%%[max_absorption_values, max_idx] = max(absorption, [], 1); % Find max absorption and its index along thickness
%%optimal_thickness_for_max_absorption = thickness_range(max_idx); % Convert index to actual thickness

% Optionally, plot the optimal thickness vs wavelength
%%figure;
%%plot(wavelength, optimal_thickness_for_max_absorption, 'LineWidth', 1);
%%xlabel('Wavelength (nm)');
%%ylabel('Optimal Thickness for Max Absorption (nm)');
%%title('Optimal Thickness vs Wavelength for Maximum Absorption');
%%grid on;



for temp_idx = 1:length(temperature_range)
    % Update waitbar less frequently to reduce overhead
    if mod(temp_idx, 5) == 0 || temp_idx == length(temperature_range)
        waitbar(temp_idx/length(temperature_range), thermal_waitbar, ...
            sprintf('Computing thermal properties: Temperature %.1f K (%.1f%%)', ...
            temperature_range(temp_idx), 100*temp_idx/length(temperature_range)));
    end
    
    T = temperature_range(temp_idx);
    
    % Pre-compute blackbody radiation for this temperature across all wavelengths
    L_bb = (2 * pi * planck_const * speed_light^2) ./ ...
           (wavelength_m.^5 .* (exp((planck_const * speed_light) ./ (wavelength_m * boltzmann_const * T)) - 1));
    
    for t_idx = 1:length(thickness_range)
        % FIXED: Calculate emissive power with explicit summation
        total_power_E = 0;
        total_power_B = 0;
        
        for w_idx = 1:length(wavelength_m)
            A = absorption(t_idx, w_idx);
            E = A * L_bb(w_idx);
            
            total_power_E = total_power_E + E * dlambda(w_idx);
            total_power_B = total_power_B + L_bb(w_idx) * dlambda(w_idx);
        end
        
        total_emissive_power_structure(t_idx, temp_idx) = total_power_E;
        total_emissive_power_structure_BB(t_idx, temp_idx) = total_power_B;
        Loss_Rate(t_idx, temp_idx) = 100 * (total_power_E / total_power_B);
    end
end

% Close the waitbar
close(thermal_waitbar);

% Display completion message
disp('Calculations complete! Generating plots...');

% % Transmission
% figure;
% mesh(wavelength, thickness_range, transmission);
% set(gca, 'YDir', 'normal');
% xlabel('Wavelength (nm)');
% ylabel('Thickness (nm)');
% title('Transmission vs Wavelength and Thickness');
% colorbar;
% clim([0 1]);
% view(2);

% % Absorption
% figure;
% mesh(wavelength, thickness_range, absorption);
% set(gca, 'YDir', 'normal');
% xlabel('Wavelength (nm)');
% ylabel('Thickness (nm)');
% title('Absorption vs Wavelength and Thickness');
% colorbar;
% caxis([0 1]);
% view(2);

% % Reflection
% figure;
% mesh(wavelength, thickness_range, reflection);
% set(gca, 'YDir', 'normal');
% xlabel('Wavelength (nm)');
% ylabel('Thickness (nm)');
% title('Reflection vs Wavelength and Thickness');
% colorbar;
% caxis([0 1]);
% view(2);


% % Emissive Power
% figure;
% mesh(temperature_range, thickness_range, total_emissive_power_structure);
% set(gca, 'YDir', 'normal');
% xlabel('Temperature (K)');
% ylabel('Thickness (nm)');
% zlabel('Total Emissive Power (W/m^2)');
% title('Total Emissive Power vs Thickness and Temperature');
% colorbar;
% view(2);


% % Loss Rate
% figure;
% mesh(temperature_range, thickness_range, Loss_Rate);
% set(gca, 'YDir', 'normal');
% xlabel('Temperature (K)');
% ylabel('Thickness (nm)');
% zlabel('Loss Rate (%)');
% title('Loss Rate vs Thickness and Temperature');
% colorbar;
% view(2);
% 
% disp('All plots have been generated successfully.');
% 


% % Plot Absorption vs Wavelength and Thickness using surf (for better control)
% figure;
% surf(wavelength, thickness_range, absorption, 'EdgeColor', 'none'); % Use surf for smooth surface
% set(gca, 'YDir', 'normal');
% xlabel('Wavelength (nm)');
% ylabel('Thickness (nm)');
% title('Absorption vs Wavelength and Thickness');
% colorbar;
% caxis([0 1]);
% view(2);
% hold on;

% Overlay the optimal thickness for max absorption in 3D by adding a small offset to the Z axis
%%z_offset = max(absorption(:)) * 1.05; % A small offset for visibility
%%plot3(wavelength, optimal_thickness_for_max_absorption, repmat(z_offset, size(wavelength)), 'r-', 'LineWidth', 0.5);

% Display the plot with the overlaid thickness curve
%legend('Absorption', 'Optimal Thickness for Max Absorption', 'Location', 'Best');





target_wavelengths = [450, 10600]; % in nm
absorption_values = zeros(size(target_wavelengths));
optimal_thicknesses = zeros(size(target_wavelengths));

for i = 1:length(target_wavelengths)
    [~, idx] = min(abs(wavelength - target_wavelengths(i)));
    [absorption_values(i), thickness_idx] = max(absorption(:, idx));
    optimal_thicknesses(i) = thickness_range(thickness_idx);
end

% Display results
fprintf('\n--- Optimal Thickness Results ---\n');
for i = 1:length(target_wavelengths)
    fprintf('Wavelength: %.0f nm -> Optimal Thickness: %.2f nm, Absorption: %.4f\n', ...
        target_wavelengths(i), optimal_thicknesses(i), absorption_values(i));
end





% Define target temperatures
target_temperatures = [1000, 1500, 2000, 2500]; % in Kelvin

% Print header
fprintf('\n--- Loss Rate at Target Temperatures for Optimal Thicknesses ---\n');

% Loop through each target wavelength
for i = 1:length(target_wavelengths)
    % Find optimal thickness for current wavelength
    [~, t_idx] = min(abs(thickness_range - optimal_thicknesses(i)));
    
    fprintf('Wavelength: %.0f nm -> Optimal Thickness: %.2f nm\n', ...
        target_wavelengths(i), optimal_thicknesses(i));
    
    % Loop through each temperature
    for j = 1:length(target_temperatures)
        % Find the index of the closest temperature
        [~, temp_idx] = min(abs(temperature_range - target_temperatures(j)));
        
        % Extract Loss Rate
        loss = Loss_Rate(t_idx, temp_idx);
        
        fprintf('  T = %d K -> Loss Rate = %.2f%%\n', target_temperatures(j), loss);
    end
    fprintf('\n'); % Add spacing between wavelengths
end




% Define threshold, target wavelengths, and target temperatures
absorption_threshold = 0.98;
target_wavelengths = [450, 10600]; % in nm
target_temperatures = [1000, 1500, 2000, 2500]; % in K

fprintf('\n--- Minimum Thickness with Absorption > %.2f and Corresponding Loss Rates ---\n', absorption_threshold);

for i = 1:length(target_wavelengths)
    % Find the index of the target wavelength
    [~, w_idx] = min(abs(wavelength - target_wavelengths(i)));
    
    % Get absorption vs thickness at this wavelength
    absorption_at_wavelength = absorption(:, w_idx);
    
    % Find first index where absorption > threshold
    valid_indices = find(absorption_at_wavelength >= absorption_threshold);
    
    if isempty(valid_indices)
        fprintf('Wavelength: %.0f nm -> No thickness found with absorption > %.2f\n\n', ...
            target_wavelengths(i), absorption_threshold);
    else
        % Get the minimum thickness and its index
        min_idx = valid_indices(1);
        min_thickness = thickness_range(min_idx);
        fprintf('Wavelength: %.0f nm -> Min Thickness = %.2f nm (Absorption = %.4f)\n', ...
            target_wavelengths(i), min_thickness, absorption_at_wavelength(min_idx));
        
        % Loop through each target temperature
        for j = 1:length(target_temperatures)
            [~, temp_idx] = min(abs(temperature_range - target_temperatures(j)));
            loss = Loss_Rate(min_idx, temp_idx);
            fprintf('  T = %d K -> Loss Rate = %.2f%%\n', target_temperatures(j), loss);
        end
        fprintf('\n'); % Space between each wavelength block
    end
end
