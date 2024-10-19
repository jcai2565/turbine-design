%% MATLAB Script to Calculate Power Output and Save to Files

% Parameters from your script
rpm = 1300; % RPM 
axial_induction = 1/3;
n_blades = 3;
density = 1.204; % kg/m^3, 20degC
dyn_visc = 1.825E-5; % kg/m-s, 20degC

% Constants
omega = rpm * 2 * pi / 60; % Convert RPM to rad/s

% VARIED PARAMETERS
u1 = 3 : 0.1 : 4.2; % m/s, length 13 vector
alpha = 4 : 0.5 : 10; % degree, length 13 vector

% GEOMETRY 
radius_vector_meter = [0.0254, 0.0334, 0.0414, 0.0495, 0.0575, 0.0655, 0.0735, 0.0815, 0.0896, 0.0976, ...
                       0.1056, 0.1136, 0.1217, 0.1297, 0.1377, 0.1457, 0.1537, 0.1618, 0.1698, 0.1778];

% Loop through each u1 and alpha value
for u = 1:length(u1)
    for a = 1:length(alpha)
        
        % Preallocate pitch and chord vectors
        pitch_vector = zeros(size(radius_vector_meter));
        chord_vector = zeros(size(radius_vector_meter));
        
        % Loop over radius_vector_meter to compute pitch and chord for each radius
        for r = 1:length(radius_vector_meter)
            R = radius_vector_meter(r);
            
            % Calculate pitch
            pitch_vector(r) = calculatePitch(axial_induction, u1(u), omega, R, deg2rad(alpha(a)));
            
            % Calculate chord
            chord_vector(r) = calculateChord(axial_induction, u1(u), omega, R, n_blades, deg2rad(alpha(a)), pitch_vector(r));
        end
        
        % Prepare the data to be written to file
        output_data = [radius_vector_meter', rad2deg(pitch_vector)', chord_vector']; % r, pitch, chord
        
        % Define file name and path
        filename = sprintf('Blades/u1_%.1f_alpha_%.1f.xlsx', u1(u), alpha(a));
        
        % Add header row and write data to .xlsx file
        header = {'r', 'pitch', 'chord'};
        writecell(header, filename, 'Sheet', 1, 'Range', 'A1');
        writematrix(output_data, filename, 'Sheet', 1, 'Range', 'A2');
        
        % Display progress
        fprintf('Data for u1 = %.1f m/s, alpha = %.1f degrees saved to %s\n', u1(u), alpha(a), filename);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%

function [cl, cd] = lookupClCd(alpha)
% Looks up the Airfoiltoolsdata

airfoil_name = 'NACA4412';
% Note: We ballpark Reynold's number (NCRIT = 5, as adviced by James & PP)
% using u ~ sqrt([1200RPM * 0.0127 m]^2 + [4 m/s]^2) and
% density ~ 1.2, dyn_visc ~ 1.82E-5. 
% Gives me about ~10000 near the middle and ~20000 near the tip. So we will
% use Re = 50000 for the calculations as the closest. 
Re = 50000;

% FILENAME TO READ FOR Cl , Cd
filename = ['airfoiltools/', airfoil_name, '-', num2str(Re)];
rawdata = readmatrix(filename);

alpha_array = rawdata(:,1);
cl_array = rawdata(:,2);
cd_array = rawdata(:,3);

% Brute force lookup -- since the data increments by 0.25, we set a
% tolerance of 0.25/2 to find "closest" alpha. If [alpha] lands between two
% tabulated rows, this will go with the lesser
tol = 0.25/2;

%Init index as a negative number so we can check if reached
closest_alpha_index = -1;

for i = 1:length(alpha_array)
    if (abs(alpha_array(i) - alpha) <= tol)
        closest_alpha_index = i;
        break;
    end
end

% alpha not found check
if (closest_alpha_index == -1)
    cl = 0;
    cd = 0;
    fprintf("Input alpha: %f was not found in the file. Zeros returned for cl, cd \n", alpha);
    return
end

cl = cl_array(closest_alpha_index);
cd = cd_array(closest_alpha_index);

end

%------------------------------------------------------------------%

function pitch = calculatePitch(a, u1, omega, R, alpha)
% Mathematical function to calculate pitch angle

u2 = (1-a)*u1;
pitch = atan(u2^2 / omega / R) - alpha;

end

%------------------------------------------------------------------%

function chord_length = calculateChord(a, u1, omega, R, n_blades, alpha, pitch)
%Mathematical function to calculate chord length

u2 = (1-a)*u1;
[cl, cd] = lookupClCd(alpha);

numerator = u1^2 * 4 * a * (1-a) * pi * R;
denominator = n_blades * (omega^2 * R^2 + (1-a)^2*u1^2) * (cl*cos(pitch + alpha) + cd*sin(pitch + alpha));

chord_length = numerator / denominator;

end

%------------------------------------------------------------------%
