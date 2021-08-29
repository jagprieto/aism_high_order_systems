clc;
clear('all');
% rng('default');
warning('off','all');
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMETERS = {};
PARAMETERS.CREATE_PDF = false;
PARAMETERS.SAMPLING_TIME = 1e-2;
PARAMETERS.PLOT_FONT_SIZE = 9.0;
PARAMETERS.CONTROL_VERSION = 0;
PARAMETERS.SYSTEM_ORDER = 3;
PARAMETERS.DISTURBANCE_TYPE = 1; % 0->No disturbance; 1->Sinusoidal disturbances ;  2->Train of sinusoidal multiple frequencies disturbance 
PARAMETERS.TOTAL_TIME = PARAMETERS.SYSTEM_ORDER+3;%ceil(0.5*PARAMETERS.SYSTEM_ORDER+8);%*PARAMETERS.SYSTEM_ORDER;

% Initial state
PARAMETERS.initial_state = zeros(PARAMETERS.SYSTEM_ORDER,1);
PARAMETERS.initial_state(1) = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors =['r','b','k','c','y','g','m'];
create_fig = true;
for step=1:7
    PARAMETERS.AISM.EPSILON = 1e-2*PARAMETERS.SAMPLING_TIME;
    PARAMETERS.AISM.color = colors(step);
    PARAMETERS.AISM.T_S = 0.5*(PARAMETERS.SYSTEM_ORDER);
    PARAMETERS.AISM.ALFA_MIN = 0.75*(PARAMETERS.SYSTEM_ORDER+1)/(PARAMETERS.AISM.T_S);
    PARAMETERS.AISM.ALFA_INITIAL = 2*PARAMETERS.AISM.ALFA_MIN;
    PARAMETERS.AISM.ALFA_MAX = 2*PARAMETERS.AISM.ALFA_INITIAL;    
    PARAMETERS.AISM.LAMBDA_INITIAL = 2*PARAMETERS.AISM.ALFA_INITIAL;  
    PARAMETERS.AISM.MU_MIN = min(PARAMETERS.SAMPLING_TIME, PARAMETERS.AISM.EPSILON*(PARAMETERS.AISM.ALFA_MAX^(PARAMETERS.SYSTEM_ORDER-1))); 
    PARAMETERS.AISM.DELTA = 0.5;    
    PARAMETERS.AISM.BETA_GAIN = 2;
    PARAMETERS.AISM.ETA_GAIN = 0.1;
    PARAMETERS.OMEGA_C_MAX = 2/PARAMETERS.AISM.ALFA_MIN;       
    PARAMETERS.AISM.KAPPA_MAX = ((PARAMETERS.AISM.MU_MIN*pi/PARAMETERS.SAMPLING_TIME)^2 + PARAMETERS.OMEGA_C_MAX)/(PARAMETERS.AISM.MU_MIN^(PARAMETERS.AISM.DELTA));    

    for k = 1:PARAMETERS.SYSTEM_ORDER
        %PARAMETERS.initial_state(k) = 2*k*2*(0.5-rand(1));
        PARAMETERS.initial_state(k) = 2*(2^k)*(0.5-rand(1));
    end
%      PARAMETERS.initial_state = 1.0e+02 *[  
%          -1.259681085792156
%   -0.011958411152515
%   -0.018747286695678
%    0.071786898222021
%    0.136676704288328
%    0.263342425808179
%   -0.381889106065859
%   -1.134100837300999
%   -0.940623730730794
%     ];
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    PARAMETERS
    ASIM = PARAMETERS.AISM
    initial = PARAMETERS.initial_state
    [SIMULATION_DATA, PARAMETERS] = run_simulation(PARAMETERS);
    plot_simulation (SIMULATION_DATA, PARAMETERS, create_fig);
    
    
    SIMULATION_DATA.AISM.system_state_history(:,k)
    
%     figure(1);
%     subplot(PARAMETERS.SYSTEM_ORDER,1,1);
%     if create_fig
%         subplot_axes = axes('Position',[.55 .77 .3 .1]);
%         box on;
%     end
%     data_size = size(SIMULATION_DATA.time,1) - 1;
%     detail_size = 200;
%     from_to_detail = data_size-detail_size:data_size;    
%     plot(subplot_axes, SIMULATION_DATA.time(from_to_detail,1), SIMULATION_DATA.AISM.system_state_history(from_to_detail,1),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
%     if create_fig
%         grid on;
%         hold on;
%         xlim([PARAMETERS.TOTAL_TIME-detail_size*PARAMETERS.SAMPLING_TIME, PARAMETERS.TOTAL_TIME]);
%     end
    create_fig = false;
%     break;
end

if PARAMETERS.CREATE_PDF
   figure(1);
   export_fig(strcat('../MANUSCRIPT/GRAPHICS/states_order_', num2str(PARAMETERS.SYSTEM_ORDER), '_disturbance_', num2str(PARAMETERS.DISTURBANCE_TYPE), '.pdf'), '-transparent', '-nocrop');
    if PARAMETERS.SYSTEM_ORDER == 3
           figure(2);
           export_fig(strcat('../MANUSCRIPT/GRAPHICS/control_order_', num2str(PARAMETERS.SYSTEM_ORDER), '_disturbance_', num2str(PARAMETERS.DISTURBANCE_TYPE), '.pdf'), '-transparent', '-nocrop');
           figure(3);
           export_fig(strcat('../MANUSCRIPT/GRAPHICS/parameters_order_', num2str(PARAMETERS.SYSTEM_ORDER), '_disturbance_', num2str(PARAMETERS.DISTURBANCE_TYPE), '.pdf'), '-transparent', '-nocrop');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation
function [SIMULATION_DATA, PARAMETERS] = run_simulation(PARAMETERS)
%     disp('...................................');
    % Simulation time
    simulation_time = 0:PARAMETERS.SAMPLING_TIME:PARAMETERS.TOTAL_TIME-PARAMETERS.SAMPLING_TIME;
    PARAMETERS.SIMULATION_STEPS = size(simulation_time, 2) + 1;
    
    % Prepare simulation data
    SIMULATION_DATA = {};
    SIMULATION_DATA.time = zeros(PARAMETERS.SIMULATION_STEPS,1);
                  
    % Adaptive integral sliding mode
    [SIMULATION_DATA.AISM, PARAMETERS] = AISM.create_simulation_data(PARAMETERS);
  
    % Nested sliding mode
%     SIMULATION_DATA.NASM = NASM.create_simulation_data(PARAMETERS);
        
    % Run simulation
    simulation_time = 0.0;
    for simulation_step = 1:PARAMETERS.SIMULATION_STEPS        
        SIMULATION_DATA.time(simulation_step) = simulation_time;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Adaptive integral sliding mode control for
        %%%% high-order uncertain nonlinear system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %------------- Disturbance -------------%  
        SIMULATION_DATA.AISM.disturbance = disturbance(simulation_time, SIMULATION_DATA.AISM.system_state, PARAMETERS);
            
        %------------- System functions -------------% 
        [aism_f,aism_b] = system_f_b(SIMULATION_DATA.AISM.system_state, PARAMETERS);
        
        %------------- Control -------------% 
        [SIMULATION_DATA.AISM.control, SIMULATION_DATA.AISM.control_state] = AISM.get_control(SIMULATION_DATA.AISM.system_state, SIMULATION_DATA.AISM.control_state, aism_f,aism_b, PARAMETERS);
                
        %------------- Save data -------------% 
        SIMULATION_DATA.AISM.system_state_history(simulation_step, :) = SIMULATION_DATA.AISM.system_state;
        SIMULATION_DATA.AISM.control_history(simulation_step, 1) = SIMULATION_DATA.AISM.control;
        SIMULATION_DATA.AISM.disturbance_history(simulation_step, :) = SIMULATION_DATA.AISM.disturbance;
        SIMULATION_DATA.AISM.control_state_history(simulation_step, :) = SIMULATION_DATA.AISM.control_state;
        
        %------------- System dynamics -------------% 
        SIMULATION_DATA.AISM.system_state = system_dynamics(SIMULATION_DATA.AISM.system_state, SIMULATION_DATA.AISM.control, SIMULATION_DATA.AISM.disturbance, PARAMETERS);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%% Adaptive integral sliding mode control for
%         %%%% high-order uncertain nonlinear system
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %------------- Disturbance -------------%  
%         SIMULATION_DATA.NASM.disturbance = disturbance(simulation_time, SIMULATION_DATA.NASM.system_state, PARAMETERS);
%         
%         %------------- System functions -------------% 
%         [nasm_f,nasm_b] = system_f_b(SIMULATION_DATA.NASM.system_state, PARAMETERS);
%         
%         %------------- Control -------------% 
%         [SIMULATION_DATA.NASM.control, SIMULATION_DATA.NASM.control_state] = NASM.get_control(SIMULATION_DATA.NASM.system_state, SIMULATION_DATA.NASM.control_state, nasm_f,nasm_b, PARAMETERS);
%                 
%         %------------- Save data -------------% 
%         SIMULATION_DATA.NASM.system_state_history(simulation_step, :) = SIMULATION_DATA.NASM.system_state;
%         SIMULATION_DATA.NASM.control_history(simulation_step, 1) = SIMULATION_DATA.NASM.control;
%         SIMULATION_DATA.NASM.disturbance_history(simulation_step, :) = SIMULATION_DATA.NASM.disturbance;
%         SIMULATION_DATA.NASM.control_state_history(simulation_step, :) = SIMULATION_DATA.NASM.control_state;
%         
%         %------------- System dynamics -------------% 
%         SIMULATION_DATA.NASM.system_state = system_dynamics(SIMULATION_DATA.NASM.system_state, SIMULATION_DATA.NASM.control, SIMULATION_DATA.NASM.disturbance, PARAMETERS);
%        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                
        % Update time
        simulation_time = simulation_time + PARAMETERS.SAMPLING_TIME;
    end
end

% Disturbance
function [disturbance] = disturbance(time, system_state, PARAMETERS)
    disturbance = zeros(PARAMETERS.SYSTEM_ORDER,1);
    if PARAMETERS.DISTURBANCE_TYPE == 1
        disturbance(PARAMETERS.SYSTEM_ORDER) = 0.1*sin(20.0*time);
    end
end

% System f and b
function [f, b] = system_f_b(system_state, PARAMETERS)
    f = zeros(PARAMETERS.SYSTEM_ORDER,1);
    f(PARAMETERS.SYSTEM_ORDER) = (system_state(2))^3.0;
    b = 1.0;
end

% System dynamics
function [system_state_new] = system_dynamics(system_state, control, disturbance, PARAMETERS)
    [f, b] = system_f_b(system_state, PARAMETERS);
    system_state_new = zeros(PARAMETERS.SYSTEM_ORDER,1);
    for k = 1:PARAMETERS.SYSTEM_ORDER-1
        system_state_new(k) = system_state(k) + (system_state(k+1) + f(k) + disturbance(k))*PARAMETERS.SAMPLING_TIME;
    end    
    k = PARAMETERS.SYSTEM_ORDER;
    system_state_new(k) = system_state(k) + (f(k) + b*control + disturbance(k))*PARAMETERS.SAMPLING_TIME;
end

% Plot simulation 
function plot_simulation(SIMULATION_DATA, PARAMETERS, create)
     
    fig1 = figure(1);
    if create
        clf(fig1);
    end
    for k = 1:PARAMETERS.SYSTEM_ORDER
        subplot(PARAMETERS.SYSTEM_ORDER,1,k);
        plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.system_state_history(:,k),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
        if create
            grid on;
            hold on;
            ylabel(strcat('$x_{', num2str(k),'}(t)$'), 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
            xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
            title(strcat('State $x_', num2str(k),'(t)$'), 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
%             xlim([PARAMETERS.TOTAL_TIME-0.7, PARAMETERS.TOTAL_TIME]);
            xlim([0.0, PARAMETERS.TOTAL_TIME]);
        end
    end
     
    fig2 = figure(2);
    if create
        clf(fig2);
    end
    plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_history,'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
    if create
        grid on;
        hold on;
        ylabel('u(t)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE, 'Interpreter','latex');
        xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        title('Control effort', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlim([0.0, PARAMETERS.TOTAL_TIME]);
    end
    
    % int_s_n, lambda, alfa, z, s, kappa, mu
    fig3 = figure(3);
    if create
        clf(fig3);
    end
    subplot(4,1,1);
    plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_state_history(:,2),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
    if create
        grid on;
        hold on;
        ylabel('$\lambda$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        title('Parameter $\lambda$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlim([0.0, PARAMETERS.TOTAL_TIME]);
    end
    subplot(4,1,2);
    plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_state_history(:,3),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
    if create
        grid on;
        hold on;
        ylabel('$\alpha$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        title('Parameter $\alpha$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlim([0.0, PARAMETERS.TOTAL_TIME]);
    end
    subplot(4,1,3);
    plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_state_history(:,6),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
    if create
            grid on;
        hold on;
        ylabel('$\kappa$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        title('Parameter $\kappa$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlim([0.0, PARAMETERS.TOTAL_TIME]);
    end
    subplot(4,1,4);
    plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_state_history(:,7),'-', 'Color', PARAMETERS.AISM.color, 'LineWidth',1.0);
    if create
        grid on;
        hold on;
    %     plot(SIMULATION_DATA.time, SIMULATION_DATA.AISM.control_state_history(:,4),'-', 'Color', 'r', 'LineWidth',1.0);
        ylabel('$\mu$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlabel('Time (s)', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        title('Parameter $\mu$', 'FontSize', PARAMETERS.PLOT_FONT_SIZE,'Interpreter','latex');
        xlim([0.0, PARAMETERS.TOTAL_TIME]);
    end
    
end