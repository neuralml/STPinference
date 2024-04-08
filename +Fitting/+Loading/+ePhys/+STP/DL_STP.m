%Class for loading Short-term Plasticity data
%

classdef DL_STP < Fitting.Loading.DataLoading
    
    properties (Constant)
        %% Possible Condition
        NONE_COND_NOTNORM = 0;
        
        %% Format
        XLS = 1;
        
        %% Noise Format
        CVS = 1;
        STDS = 2;        
    end    
    
    properties
        freqs = [30];
        rec_pulses = 0;
        noise_format = 2;
        
        %% XLS Columns
        INDEX_FILE = 1;
        INDEX_QUALITY = 2;
        INDEX_CHECK = 3;
        INDEX_TYPE = 4;
        INDEX_CONDITION = 5;
        INDEX_DRUG = 6;
        INDEX_AREA = 7;
        INDEX_FROM = 8;
        
        %% NUM col
        INDEX_AGE = 16;
        INDEX_CV = 17;
    end
    
    methods
        
        %
        %Contructor
        %
        function dl = DL_STP(condition)
            disp('-- Created DL_STP --');
            
            if nargin == 0
                input = {};
            else
                input = {condition};
            end
            
            dl = dl@Fitting.Loading.DataLoading(input); % call superclass constructor
            dl.format = dl.XLS;
            dl.noise_format = dl.STDS;
        end
        
        function [d CVs files STDs list_freqs] = load(obj, path)
            STDs = [];
            CVs = [];
            
            %% 1. Load file
            [num txt] = loadXLS(path, obj);
            
            %Remove first row
            if(isnan(num(1,1)))
                num(1,:) = [];    
            end    
            %cols = txt(1,:); %Stores the columns
            
            txt(1,:) = [];
            num(:,2) = [];
            
            %Remove non-desired quality
            for a=1:length(obj.remove_quality)
                num(strcmpi(txt(:,obj.INDEX_QUALITY), obj.remove_quality{a}), :) = []; %Remove from num table
                txt(strcmpi(txt(:,obj.INDEX_QUALITY), obj.remove_quality{a}), :) = []; %Remove from chars table
            end

            %Remove non-desired drugs (all the drugs are desirable, though)
            for i=1:length(obj.remove_drugs)
                num(strcmpi(txt(:,obj.INDEX_DRUG), obj.remove_drugs(i)),:) = [];
                txt(strcmpi(txt(:,obj.INDEX_DRUG), obj.remove_drugs(i)),:) = [];
            end

            %/num(strcmpi(txt(:,6),remove_areas(1)),:) = [];
            
            %Remove entries without min number of pulses
            if(obj.min_pulses)
                txt(isnan(num(:, 1+obj.min_pulses)),:) = [];
                num(isnan(num(:, 1+obj.min_pulses)),:) = [];
            end
            
            %Remove entries with a diff age
            if(obj.age>=0)
                num(~ismember(num(:,obj.INDEX_AGE),obj.age),:) = [];
            end    
            
            %Find out how many diff frequencies
            %ufreqs = sort(unique(num(:,1)));
            %n_freq = length(freqs);
            
            list_freqs = num((ismember(num(:,1), obj.freqs)), 1);
            %keyboard
            
            switch obj.condition
                case obj.NONE_COND_NOTNORM
                    d = num((ismember(num(:,1), obj.freqs)), 2:(1+obj.min_pulses));
                    
                    if(obj.noise_format==obj.CVS) %Deals with stored CVs or STDs
                        CVs = num((ismember(num(:,1), obj.freqs)), obj.INDEX_CV:obj.INDEX_CV+obj.min_pulses-1);
                    elseif(obj.noise_format==obj.STDS)
                        STDs = num((ismember(num(:,1), obj.freqs)), obj.INDEX_CV:obj.INDEX_CV+obj.min_pulses-1);
                        CVs = zeros(size(STDs));
                    end    
                    
                    files = txt((ismember(num(:,1), obj.freqs)), obj.INDEX_FILE);
                    for i=1:size(d,1)  
                        if(obj.noise_format==obj.STDS)
                            CVs(i,:) = STDs(i,:)./d(i,1); %Normalizes to the first peak
                        end                   
                    end
            end
            files
            d
        end
       
        
        
        function [spikes stimes] = setInput(dl, stime)
            
            %Defines the input for the model being optimized
            if(dl.rec_pulses)
                [spikes stimes] = inhreg(stime, dl.freqs(1), 0);
                [spikes2 stimes2] = recovery_pulses(dl.rec_pulses);
                spikes = [spikes spikes2];
                stimes = [stimes stimes2];
            else
                [spikes stimes] = inhreg(stime, dl.freqs(1), 0);
            end        
        end  
        
    end
    
end

