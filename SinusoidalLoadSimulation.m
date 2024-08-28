% CC BY-SA (Attribution Sharealike) Carsen Villegas XXX 2024
% Written by Carsen (carsen3117@gmail.com) for Cynthia
% (cynthia.furse@utah.edu)

% Using FDTD Simulation to replicate a SHORT/MATHCED/OPEN/RESISTIVE load

% TODO:
% envelope : maximum absolute value of a standing wave. OVERLAY A DASHED LINE ON THE "TOTAL"
% GRAPH. DASHED LINE ENVELOPE OVERLAYED OVER THE TOTAL STANDING WAVE OR
% SOMETHING. MAKE A SEPARATE ENVELOPE FIGURE TO LABLE LATER. SHORT OPEN
% MATCH AND RANDOM RESISTOR OR TWO (4 FIGURES).  the envelope allows you to label them. Use powerpoint with
% each "video" in a powerpoint to embed in the powerpoint. Name of each

% file associated with the simulation.  Two separate simulations : (one for steps and pulses) one for
% (sine wave). Can use same programming for step and pulse. 

clear all;
close all;

global load;

% SETUP------------------
% Physical Constants
eps_0 = 8.85418782E-12;  %  Permittivity of free space (F/m).
mu_0 = 4.0*pi*1E-7;      %  Permeability of free space (H/m).

% Simulation parameters
dz = 0.001;                % Grid resolution (m)
kmax  = 1000;              % Number of grid points to model 1 m of transmission line
nmax  = 2500;              % Number of time steps

% INCIDENT WAVE DATA SET
v = zeros(1, kmax);
i = zeros(1, kmax - 1);

% REFLECTED WAVE DATA SET
v_reflected = zeros(1, kmax);
i_reflected = zeros(1, kmax - 1);

% ENVELOPE DATA
envelope_data = zeros(1,kmax);


% Source parameters
Vg = 1.0;                   % Amplitude (V)
freq = 1e9;                 % Frequency for the source 1 GHz

% Transmission line parameters for RG-58
b   = (3.51e-3)/2;       % Outer conductor radius (m)
a   = (0.91e-3)/2;       % Inner conductor radius (m)

% Conductor information (copper)
sigma_cond = 5.797*1E7; % Conductivity (S/m) of copper
epsr_cond = 1.0;        % Relative permittivity of copper
mu_cond = mu_0;         % Permeability of copper

% Insulator (dielectric) information (polyethylene)
sigma_ins = 0.0001;     % Conductivity (S/m)
epsr_ins = 2.25;        % Relative permittivity
mu_ins = mu_0;          % Permeability of the insulator

% Solving for RLGC constants of coaxial cable:
Rs = sqrt(pi * freq * mu_cond / sigma_cond);
Rprime = (Rs / (2 * pi)) * ((1 / a) + (1 / b));
Lprime = (mu_ins / (2 * pi)) * log(b / a);
Gprime = (2 * pi * sigma_ins) / log(b / a);
Cprime = (2 * pi * epsr_ins * eps_0) / log(b / a);

% Angular frequency
w = 2 * pi * freq;

% Prop. constant
gamma = sqrt((Rprime + 1i * w * Lprime) * (Gprime + 1i * w * Cprime));

% Wave propagation speed
up = w / imag(gamma);

% Time step increment
dt = dz / up;

% FDTD coefficients:
c1 = -(2 * dt) / (dt * dz * Rprime + 2 * dz * Lprime);
c2 = (2 * Lprime - dt * Rprime) / (2 * Lprime + dt * Rprime);
c3 = -(2 * dt) / (dt * dz * Gprime + 2 * dz * Cprime);
c4 = (2 * Cprime - dt * Gprime) / (2 * Cprime + dt * Gprime);
%--------------------------------------------------------------

% PLOTTING

% Plot settings
Vmax_plot = 2;            % Max y-scale (V) for the voltage plot
Imax_plot = 40;           % Max y-scale (mA) for the current plot
T_fr = 2;                 % Time steps per plot update

% Set the x-axis (m) in the voltage plot
z_voltage  = (0:kmax-1) * dz;
z_current  = (0:kmax-2) * dz + dz / 2;

% Initialize a figure with three subplots
fig = figure(1);
clf;
set(fig, 'units', 'pixels', 'position', [100, 100, 640, 480]);
set(fig, 'color', [1 1 1] / 1.5);

% INCIDENT WAVE
subplot(3, 1, 1);
lineVIncident = plot(z_voltage, v, 'b', 'linewidth', 2);  % Blue for source wave
axis([min(z_voltage) max(z_voltage) -Vmax_plot Vmax_plot]);
grid on;
ylabel('Incident (V)');
xlabel('Position (m)');
set(gca, 'FontSize', 14);

% REFLECTED WAVE
subplot(3, 1, 2);
lineVReflected = plot(z_voltage(1:end-1), v_reflected(1:end-1), 'r', 'linewidth', 2);  % Red for reflected wave
axis([min(z_voltage) max(z_voltage) -Vmax_plot Vmax_plot]);
grid on;
ylabel('Reflected (V)');
xlabel('Position (m)');
set(gca, 'FontSize', 14);

% TOTAL WAVE // ENVELOPE
subplot(3, 1, 3);
lineV = plot(z_voltage(1:end-1), v(1:end-1) + v_reflected(1:end-1), 'k', 'linewidth', 2);  % Black for total wave
hold on
%Envelope below
Envelope = plot(z_voltage(1:end-1), zeros(1,kmax-1), 'k', 'Linestyle', ':', 'linewidth', 2, 'visible', 'off');
axis([min(z_voltage) max(z_voltage) -Vmax_plot Vmax_plot]);
grid on;
ylabel('Total (V)');
xlabel('Position (m)');
set(gca, 'FontSize', 14);
hold off

steady_state = false;

% EXTRA FIGURE STUFF---------------------INTERACTIVITY-------------
load = 'SHORT';
sgtitle(strcat("Transmission for ", load, " load (1GHz)"));
global m
m = 1; % counter for envelope

% Add buttons for interactivity - PLAY, PAUSE, RESTART, AND SHOW ENVELOPE
uicontrol('Style', 'pushbutton', 'String', 'Pause', 'Position', [10 10 50 20], 'Callback', 'pause_sim = true;');
uicontrol('Style', 'pushbutton', 'String', 'Resume', 'Position', [70 10 50 20], 'Callback', 'pause_sim = false;');
uicontrol('Style', 'pushbutton', 'String', 'Restart', 'Position', [130 10 50 20], 'Callback', 'restart_sim = true;');
uicontrol('Style', 'pushbutton', 'String', 'Toggle Envelope', ...
    'Position', [450 10 100 20], ...
    'Callback', @(src,event) toggleVisibility(Envelope));


% BUTTONS FOR CHANGING CONDITIONS
% Create the popup menu for boundary conditions
uicontrol('Style', 'popupmenu', 'String', {'SHORT', 'OPEN', 'MATCH', 'RESISTIVE (27 Ω)'}, ...
    'Position', [190, 10, 80, 20], ...
    'Callback', @popupMenuCallback);



% Ensure RL for "Incident" is matched
Z0 = sqrt((Rprime + 1i * w * Lprime) / (Gprime + 1i * w * Cprime));

pause_sim = false;
restart_sim = false;
show_envelope = false; % TODO
global reset_envelope
reset_envelope = false;

n = 1;
% Start time-stepping---------------------
while true
    % Simulation control:
    if restart_sim
        v = zeros(1, kmax);
        i = zeros(1, kmax - 1);
        v_reflected = zeros(1, kmax);
        i_reflected = zeros(1, kmax - 1);
        restart_sim = false;
        envelope_data = zeros(1,kmax);
        n = 1;
        m = 1;
    end

    if reset_envelope
        envelope_data = zeros(1,kmax);
        set(Envelope, 'ydata', envelope_data, 'visible', 'off');
        reset_envelope = false;
        m = 1000;
    end

    while pause_sim
        pause(0.1);
        continue;
    end

    % TODO : FIX frequency here actually being something good. 
    v(1) = sin(n*dt * 1e9 * 2*pi);

    % Right End Boundary Condition (MATCHED for "INCIDENT")
    v(kmax) = abs(Z0) * i(kmax - 1);
    i(kmax - 1) = v(kmax) / abs(Z0); 

    % Left End boundary Condition (MATCHED for "REFLECTED")
    v_reflected(1) = v_reflected(2) - abs(Z0) * (i_reflected(1) - i_reflected(2));


%--------------------

    % Simulation updating:
    % Current updating loop
    for k = 1:kmax - 1
        i(k) = c1 * (v(k + 1) - v(k)) + c2 * i(k);
        i_reflected(k) = c1 * (v_reflected(k + 1) - v_reflected(k)) + c2 * i_reflected(k);
    end

    % Voltage updating loop
    for k = 2:kmax - 1
        v(k) = c3 * (i(k) - i(k - 1)) + c4 * v(k);
        v_reflected(k) = c3 * (i_reflected(k) - i_reflected(k - 1)) + c4 * v_reflected(k);
        v_reflected(kmax) = vReflectedHandling(load, v(kmax)); % update boundary according to boundary conditions
    end

    % Envelope updating loop
    if m > 2000
        show_envelope = true;
        for k = 1:kmax
            if abs(v(k) + v_reflected(k)) > envelope_data(k)
                envelope_data(k) = abs(v(k) + v_reflected(k));
            end
        end
    end

    

    % n = 2000 for incident and reflected wave to reach steady state
    % n = 1000 for a wave to reach the other end



    % Only update the plot every 2 time steps so it will plot faster
    if mod(n, T_fr) == 0
        set(lineVIncident, 'ydata', v);
        set(lineVReflected, 'ydata', v_reflected(1:end-1));
        v_total = v(1:end-1) + v_reflected(1:end-1);
        set(lineV, 'ydata', v_total);
        if m == 2000
            set(Envelope, 'visible', 'on')
        end
        set(Envelope, 'ydata', envelope_data(1:end-1));
        drawnow;
    end
    n = n + 1;
    m = m + 1;
end



% FUNCTIONS-----------------------------------------

function popupMenuCallback(src, ~)
    % Access the global variable
    global load;
    global reset_envelope
    reset_envelope = true;
    
    % Get the selected index
    selectedIndex = src.Value;
    
    % Define the options corresponding to the index
    options = {'SHORT', 'OPEN', 'MATCH', 'RESISTIVE (27 \Omega)'};
    
    % Update the global variable based on the selection
    load = options{selectedIndex};
    
    sgtitle(strcat("Transmission for ", load, " load (1GHz)"));

    
    % Optional: Display the selected load type
    disp(['Selected load type: ', load]);
end

% This function handles different types of reflected loads and returns
% the appropriate value based on the load type.
% Switch case to handle different load types
function vReflected = vReflectedHandling(loadType, incidentVoltage)

    switch loadType
        case 'MATCH'
            % Matched Load: vReflected should be 0 (no reflection)
            vReflected = 0;

        case 'OPEN'
            % Open Load: vReflected should be 0 (open circuit)
            vReflected = 1 * incidentVoltage;

        case 'SHORT'
            % Short Load: vReflected should be -i * Z0 (short circuit reflection)
            vReflected = -1 * incidentVoltage;

        case 'RESISTIVE (27 \Omega)'
            % Resistive Load: Example resistance value
            vReflected = .33 * incidentVoltage;

        otherwise
            % Handle unknown load types
            error('Unknown load type: %s', loadType);
    end
end

function toggleVisibility(Envelope)
    global m;
    if m > 2000
        if get(Envelope, 'Visible')
            set(Envelope, 'Visible', 'off')
        else
            set(Envelope, 'Visible', 'on')
        end
    end

end



