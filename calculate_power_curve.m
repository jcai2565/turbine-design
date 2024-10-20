function [P, Cp] = calculate_power_curve(U1, rho, RPM, n_blades, n, geometry_file, airfoil_file)
% CALCULATE_POWER_CURVE Computes theoretical power curve for a small-scale wind turbine
%
% INPUTS:
%   U1             - Free stream velocity [m/s]
%   rho            - Air density [kg/m^3]
%   RPM            - Desired RPM
%   n_blades       - Number of blades
%   n              - Number of airfoil sections to discretize blade
%   geometry_file  - Path to blade geometry file (.xlsx)
%   airfoil_file   - Path to airfoil data file (.xlsx)
%
% OUTPUTS:
%   P              - Total power extracted by the turbine [W]
%   Cp             - Power coefficient [-]

    %% Read and process user-provided data
    % Read blade geometry data
    geometry = readtable(geometry_file);

    % Read CL and CD data for the airfoil
    cl_cd_data = readtable(airfoil_file, 'Sheet', 'Sheet1');

    % Blade dimensions (from geometry data)
    R_root = geometry.r(1);      % Root radius [m]
    R_tip = geometry.r(end);     % Tip radius [m]

    % Convert RPM to rad/s
    Omega = (RPM/60) * 2 * pi; 

    %% Discretize blade into sections
    % Create radial locations
    r0 = linspace(R_root, R_tip, n+1);
    dr = diff(r0);                % Section spans
    r = (r0(1:end-1) + r0(2:end)) / 2;  % Section midpoints

    % Interpolate chord and pitch for each section
    chords = interp1(geometry.r, geometry.chord, r, 'spline');
    pitches = interp1(geometry.r, geometry.pitch, r, 'spline');

    %% Initialize variables for power calculation
    dT = zeros(1, n);  % Axial force
    dQ = zeros(1, n);  % Useful torque
    aval = zeros(1, n);  % Axial induction factor
    apval = zeros(1, n); % Angular induction factor
    alphaval = zeros(1, n); % Angle of attack
    clval = zeros(1, n);   % Lift coefficient
    cdval = zeros(1, n);   % Drag coefficient

    % Minimization parameters
    options = optimset('TolX', 1e-10, 'MaxFunEvals', 10000000);
    crit = 1e-3;  % Convergence criterion

    %% Power calculation loop for each blade section
    for i = 1:n
        chord = chords(i);
        pitch = pitches(i);

        % Initial guess for a and a_prime
        if i == 1
            as0 = [1/3, 0.0];  % Initial guess for the first section
        else
            as0 = [aval(i-1), apval(i-1)];  % Use previous values
        end

        % Find optimal a and aprime using fminsearch
        [asmin, fval, exitflag] = fminsearch(@(as)getF(as(1), as(2), chord, pitch, r(i), dr(i), cl_cd_data, U1, rho, Omega, n_blades), as0, options);

        % Store values
        aval(i) = asmin(1);
        apval(i) = asmin(2);

        % Calculate angle of attack and lift/drag coefficients
        u = U1 * (1 - aval(i));
        w = Omega * r(i) * (1 + apval(i));
        Urel = sqrt(u^2 + w^2);
        phi = atan(u / w);
        alphaval(i) = rad2deg(phi) - pitch;
        clval(i) = interp1(cl_cd_data.alpha, cl_cd_data.cl, alphaval(i), 'spline');
        cdval(i) = interp1(cl_cd_data.alpha, cl_cd_data.cd, alphaval(i), 'spline');

        % Check for valid convergence
        if fval < crit && aval(i) < 1
            dT(i) = rho * U1^2 * 4 * aval(i) * (1 - aval(i)) * pi * r(i) * dr(i);
            dQ(i) = 4 * apval(i) * (1 - aval(i)) * rho * U1 * pi * Omega * r(i)^3 * dr(i);
        else
            fprintf('MODEL CONVERGENCE FAILURE at radial location r=%f, STOPPING CALCULATION\n', r(i));
            dT = 0;
            dQ = 0;
            break;
        end
    end

    %% Evaluate total power and power coefficient
    P = sum(dQ) * Omega;  % Total power
    Cp = P / (0.5 * rho * pi * R_tip^2 * U1^3);  % Power coefficient

    %% Plot and save results
    % Create 2x4 subplot for key parameters without displaying
    fig = figure('visible', 'off');
    
    subplot(2, 4, 1);
    plot(r, clval);
    title('Cl');
    xlabel('Radius [m]');
    ylabel('Cl [-]');

    subplot(2, 4, 2);
    plot(r, cdval);
    title('Cd');
    xlabel('Radius [m]');
    ylabel('Cd [-]');

    subplot(2, 4, 3);
    plot(r, alphaval);
    title('Angle of attack');
    xlabel('Radius [m]');
    ylabel('Alpha [deg]');

    subplot(2, 4, 4);
    plot(r, aval);
    title('Axial induction factor');
    xlabel('Radius [m]');
    ylabel('a [-]');

    subplot(2, 4, 5);
    plot(r, apval);
    title('Angular induction factor');
    xlabel('Radius [m]');
    ylabel('aprime [-]');

    subplot(2, 4, 6);
    plot(r, dT ./ dr);
    title('Axial force');
    xlabel('Radius [m]');
    ylabel('Axial force per unit length');

    subplot(2, 4, 7);
    plot(r, dQ ./ (r .* dr));
    title('Useful force');
    xlabel('Radius [m]');
    ylabel('Useful force per unit length');

    subplot(2, 4, 8);
    plot(r, Omega * dQ ./ dr);
    title('Power');
    xlabel('Radius [m]');
    ylabel('Power per unit length');

    % Save the plot as an image with identifier
    if ~exist('Plots', 'dir')
        mkdir('Plots');  % Create Plots directory if it doesn't exist
    end
    [~, geom_name, ~] = fileparts(geometry_file);  % Extract geometry filename
    saveas(fig, fullfile('Plots', sprintf('power_curve_plot_%s.png', geom_name)));
    close(fig);  % Close the figure after saving

end

function [F] = getF(a, aprime, chord, pitch, r, dr, cl_cd_data, U1, rho, Omega, n_blades)
    % Evaluate thrust/torque from momentum conservation
    if a < 0.4
        dTm = rho * U1^2 * 4 * a * (1 - a) * pi * r * dr;  % Low-induction model
    else
        dTm = rho * U1^2 * (0.889 - (0.0203 - (a - 0.143)^2) / 0.6427) * pi * r * dr;  % Glauert model
    end
    dQm = rho * U1 * 4 * aprime * (1 - a) * pi * Omega * r^3 * dr;

    % Velocity components
    u = U1 * (1 - a);
    w = Omega * r * (1 + aprime);
    Urel = sqrt(u^2 + w^2);
    phi = atan(u / w);
    alpha = phi - deg2rad(pitch);

    % Lift and drag coefficients
    cl = interp1(cl_cd_data.alpha, cl_cd_data.cl, rad2deg(alpha), 'spline');
    cd = interp1(cl_cd_data.alpha, cl_cd_data.cd, rad2deg(alpha), 'spline');

    % Lift and drag
    dFL = 0.5 * rho * Urel^2 * cl * chord * dr;
    dFD = 0.5 * rho * Urel^2 * cd * chord * dr;

    % Thrust and torque from airfoil theory
    dTa = n_blades * (dFL * cos(phi) + dFD * sin(phi));
    dQa = n_blades * (dFL * sin(phi) - dFD * cos(phi)) * r;

    % Cost function
    F = (dTa - dTm)^2 + (dQa - dQm)^2;
end
