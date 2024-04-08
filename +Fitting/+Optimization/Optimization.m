%% Implements super class optimization
%

classdef Optimization < handle
    
    properties
        min_error = 0.0001;
        eval; %Evaluation function
        p_on; %Parameters being optimized
        p_limits;
        online_plot = 0;
        runs = 1;
        norm2before = 0;
        cond_reuse = 0; %Starts the opt with the best from the other cond
    end
    
    methods
        
        %
        %Contructor
        %
        function opt = Optimization()
        end
        
        
        function params = run(opt, data)
            
        end 
        
        
    end
    
end

