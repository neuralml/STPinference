classdef LearningRule < matlab.mixin.Copyable
    
    properties
         %Ref to pre/post neurons
         preneuron;
         postneuron;
         winit;
         
         run_fun;
    end
    
    
    methods
        %
        %Contructor
        %
        function lr = LearningRule()
           
        end
        
        function reset(lr)
        end    
        
        function initRule(syn, time, n, winit, dt)
        end
        
        function setBounds(lr, wmin, wmax)
        end 
        
        function plot(lr, time)
        end    
      
    end
    
end

