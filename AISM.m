classdef AISM
   methods (Static)       
       function [SIMULATION_DATA, PARAMETERS_new] = create_simulation_data(PARAMETERS)
        SIMULATION_DATA = {};       
        SIMULATION_DATA.system_state_history = zeros(PARAMETERS.SIMULATION_STEPS, PARAMETERS.SYSTEM_ORDER); 
        SIMULATION_DATA.disturbance = 0.0;
        SIMULATION_DATA.disturbance_history = zeros(PARAMETERS.SIMULATION_STEPS, PARAMETERS.SYSTEM_ORDER); 
        SIMULATION_DATA.control = 0.0;
        SIMULATION_DATA.control_history = zeros(PARAMETERS.SIMULATION_STEPS, 1);     
        SIMULATION_DATA.control_state_history = zeros(PARAMETERS.SIMULATION_STEPS, 10); 
       [SIMULATION_DATA.system_state, SIMULATION_DATA.control_state, PARAMETERS_new] = AISM.initialize_control(PARAMETERS.initial_state, PARAMETERS);    
      end 
       
      function [system_state_new, control_state, PARAMETERS_new]  = initialize_control(system_state, PARAMETERS)
        PARAMETERS_new = PARAMETERS;      
        system_state_new = system_state;
        control_state = zeros(10, 1);   % int_s_n, lambda, alfa, z, s, kappa, mu
        control_state(2,1) = PARAMETERS.AISM.LAMBDA_INITIAL; % lambda
        control_state(3,1) = PARAMETERS.AISM.ALFA_INITIAL; % alfa
        [s_n, p_n] = AISM.get_s_n(system_state, control_state(3,1), PARAMETERS.SYSTEM_ORDER);
        control_state(5,1) = s_n;
      end
      
      function [p_n] = get_p_n(system_state, alfa, n)
        p_n = 0.0;
        if n > 1
            for k=1:n-1
                p_n = p_n +  nchoosek(n-1,k)*(alfa^k)*system_state(n-k);
            end
        end
      end
      
      function [s_n, p_n] = get_s_n(system_state, alfa, n)
        if n == 1
            p_n = 0.0;
        else            
            p_n = AISM.get_p_n(system_state, alfa, n);            
        end
        s_n = system_state(n) + p_n;
      end
      
      function [dot_s_n] = get_dot_s_n(system_state, alfa, dot_alfa, n)
        if n == 1
            p_n = 0.0;
            s_n_1 = 0.0;
        else            
            p_n = AISM.get_p_n(system_state, alfa, n);  
            [s_n_1, p_n_1] = AISM.get_s_n(system_state, alfa, n-1);        
        end
        dot_s_n = system_state(n+1) + p_n + (n-1)*dot_alfa*s_n_1;
      end
                  
      function [control, control_state_new]  = get_control(system_state, control_state, f, b, PARAMETERS) 
        
        % Compute Sn and Sns-1
        alfa = control_state(3,1);
        [s_n_1, p_n_1] = AISM.get_s_n(system_state, alfa, PARAMETERS.SYSTEM_ORDER-1);
        [s_n, p_n] = AISM.get_s_n(system_state, alfa, PARAMETERS.SYSTEM_ORDER);
        
        % Compute z
        int_s_n = control_state(1,1);
        lambda = control_state(2,1);
        gamma = ((lambda)^2)/4.0; 
        z = s_n + (lambda/2.0)*int_s_n; 
        
        % Compute mu        
        mu =  max(PARAMETERS.AISM.MU_MIN, PARAMETERS.AISM.ETA_GAIN*tanh(abs(z)));
        
        % Compute alfa
        dot_alfa = -(PARAMETERS.AISM.BETA_GAIN*mu/(PARAMETERS.SYSTEM_ORDER-1))*(sign(s_n)/s_n_1);
        alfa = alfa + dot_alfa*PARAMETERS.SAMPLING_TIME; 
        if alfa > PARAMETERS.AISM.ALFA_MAX
             alfa  = PARAMETERS.AISM.ALFA_MAX;
             dot_alfa = 0;
        elseif alfa < PARAMETERS.AISM.ALFA_MIN
             alfa  = PARAMETERS.AISM.ALFA_MIN;
             dot_alfa = 0;
        end       
        
        % Compute kappa and lambda
         if abs(z) < mu
            dot_lambda = 0;
            kappa = control_state(6,1);
         else
            omega_c = 2.0/alfa;
            kappa = min(PARAMETERS.AISM.KAPPA_MAX, omega_c/(mu^((1+PARAMETERS.AISM.DELTA))));       
            dot_lambda = kappa*((abs(z))^(PARAMETERS.AISM.DELTA))*sign(z)*sign(s_n);
        end
        lambda = lambda + dot_lambda*PARAMETERS.SAMPLING_TIME;
        
        % Compute control
        control = (-1.0/b)*(f(PARAMETERS.SYSTEM_ORDER) + p_n + 1*(PARAMETERS.SYSTEM_ORDER-1)*dot_alfa*s_n_1 + lambda*s_n  + gamma*int_s_n);  
        
        % Update control state
        control_state_new = control_state;        
        control_state_new(1,1) = control_state_new(1,1) + s_n*PARAMETERS.SAMPLING_TIME;  % int_s_n
        control_state_new(2,1) = lambda;  
        control_state_new(3,1) = alfa; 
        control_state_new(4,1) = z; 
        control_state_new(5,1) = s_n; 
        control_state_new(6,1) = kappa;        
        control_state_new(7,1) = mu; 
      end
   end
end
