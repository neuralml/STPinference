%Super class for DataLoading classes
%

classdef DataLoading < handle
    
    properties
        min_pulses;
        condition;
        format;
        age = -1; %Select a given age interval
        remove_quality = {'Poor'};
        remove_areas = {''};
        remove_drugs = {''};
        
        stds; %Standard deviations
    end
    
    methods
        
        %
        %Contructor
        %
        function dl = DataLoading(input)
            dl.condition = input{1};
        end
        
    end
    
end

