% Define constants
lambda = 5;             % Optimal tip speed ratio (constant)
maxR = 0.1778;          % Max of radius_vector_meter from the generation
rho = 1.204;            % Air density [kg/m^3]
n_blades = 3;           % Number of blades
n = 100;                % Number of airfoil sections
airfoil_file = 'airfoiltools/NACA4412-50000.xlsx';  % Airfoil file
output_file = 'Data/output_results.xlsx';  % Output Excel file

% Define the directory containing geometry files
geometry_dir = 'Blades/';
files = dir(fullfile(geometry_dir, 'u1_*.xlsx'));  % Get all geometry files

% Initialize the results table
results = cell(length(files) + 1, 4);
results(1, :) = {'alpha', 'u1', 'P', 'Cp'};  % Set header row

% Loop through each geometry file
for i = 1:length(files)
   % Get the filename and extract u1 and alpha from the filename
   filename = files(i).name;
   file_path = fullfile(geometry_dir, filename);
  
   % Extract u1 and alpha from the filename
   tokens = regexp(filename, 'u1_(\d+(\.\d+)?)_alpha_(\d+(\.\d+)?)', 'tokens');
   u1 = str2double(tokens{1}{1});    % Extract u1 value
   alpha = str2double(tokens{1}{2}); % Extract alpha value
  
   % Call the calculate_power_curve function
   fprintf('Progress: %d/%d...\n', i, length(files));
   [P, Cp] = calculate_power_curve(u1, rho, (omega_calculation(lambda, maxR, u1)*60/2/pi), n_blades, n, file_path, airfoil_file);
  
   % Store the results in the cell array
   results{i + 1, 1} = alpha;
   results{i + 1, 2} = u1;
   results{i + 1, 3} = P;
   results{i + 1, 4} = Cp;
end

% Write the results to an Excel file
writecell(results, output_file);

% Display completion message
fprintf('Results saved to %s\n', output_file);

function omega = omega_calculation(lambda, R, u)
%Mathematical function to calculate optimal rotational rate to achieve tip
% speed ratio = [lambda]
   omega = lambda*u/R;
end

