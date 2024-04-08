%Class with analysis functions
%

classdef Analyser < handle
    
    properties
        sign_level = 0.05;
        fontsize = 11;
        dl;
    end
    
    methods
        
        %
        %Contructor
        %
        function an = Analyser()
            
        end
        
        function results = distPlot(an, names, p_limits, varargin)
            points = 1000;
            
            if (nargin>1)
                dens = varargin{1};
                best1 = varargin{2};
            end    
            
            figure();
            %slopes = zeros(size(params1,1),1);
            results1 = cell(size(dens,1), length(names)+1);
            for i=1:length(names)                
                if(length(names)>1)
                    subplot(round(length(names)/min(length(names),2)), min(length(names),2), i);
                end
                c = 1;
                y_mean = zeros(1, points);
                xbins = linspace(p_limits(i,1), p_limits(i,3), points);
                for j=1:size(dens,2)
                    
                    p = marginal(dens{j}, i);
                    y = evaluate(p,xbins);
                    
                    plot(xbins, y, 'Color', [0.9 0.5 0.5], 'LineWidth', 1.2);
                    hold on;                    
                    ylabel('Posterior Probability', 'Fontsize', 11);
                    y_mean = y_mean + y;
                    
                    scatter(best1(j,i), evaluate(p, best1(j,i)), 45, [0.9 0.5 0.5],'filled');
                    c=c+1;
                end
                %if(size(f_mean2,1)>1)
                    plot(xbins, y_mean./size(dens,2), 'Color', [1 0 0], 'LineWidth', 2);
                %end    
                xlim([p_limits(i,1) p_limits(i,3)]);
                
                n = names{i};
                
                xlabel(n, 'Fontsize', 11);
                box off;
                hold off;
            end
            set(gcf, 'Color', 'w');
            results = results1;
        end 
        
        
        function [results ppr rsq rsq_tot adj_rsq adj_rsq_tot pred peaks] = compareResults(an, data, names, params, models, p_on, norm2before)
            %Analysis if there is any statistical significant difference
            %between the data and the simulation
            
            %pre_peaks = [];
            results = zeros(size(params,1),2);
            ppr = zeros(size(params,1),2);
            pred = zeros(size(data));
            figure;
            peaks_bef = [];
            for i=1:size(params,1)
                
                %% Peaks fitting, TODO: EPSP fitting?                
                
                if(norm2before==3) %Before/after at the same time
                    models(i).setParams(params(i,:), p_on);
                    models(i).reset();
                    [peaks_before peaks_after] = models(i).run_fun(); %Runs model
                    peaks = [peaks_before./peaks_before(1) peaks_after./peaks_before(1)];
                    peaks_bef = peaks;
                else    
                    models(i).setParams(params(i,:), p_on);
                    models(i).reset();
                    [peaks] = models(i).run_fun(); %Runs model
                    peaks_bef2 = peaks;
                    
                    if(norm2before && mod(i,2) == 0)
                        %peaks = peaks.*(dot(data(i-1,:),peaks./peaks_bef(1))/norm(peaks./peaks_bef(1))^2); %Scale by A_MAP before
                        peaks = peaks.*sum(data(i-1,:).*(peaks_bef)./(an.dl.stds(i,:).^2))/sum(((peaks_bef).^2)./(an.dl.stds(i,:).^2)); %Scale by A_MAP
                        %peaks = peaks./pre_peaks(1);
                    else
                        %pre_peaks = peaks;
                        %peaks = peaks./peaks(1);                    
                        peaks = peaks.*(dot(data(i,:),peaks)/norm(peaks)^2); %Scale by A_MAP
                    end
                    peaks_bef = peaks_bef2;
                end
                ppr(i,:) = [(data(i,2)-data(i,1))/data(i,1) (peaks(2)-peaks(1))/peaks(1)];
                
                pred(i,:) = peaks;
                [h,p,ci] = ttest(data(i,:)./data(i,1), peaks./peaks(1), 0.8,'both');
                results(i,:) = [h p];
                
                subplot(ceil(size(params,1)/2), 2, i);
                
                scatter(1:length(peaks), peaks*1e3,'r','filled');
                hold on;
                scatter(1:length(peaks), data(i,:)*1e3,'g','filled');
                xlabel('Pulse number', 'fontsize', an.fontsize-2);
                ylabel('PSP (mV)', 'fontsize', an.fontsize-2);
                
                title(strcat(names{i}, ' f=', num2str((models(i).freq)), 'Hz'), 'fontsize', an.fontsize-1);
                if(i==1)
                    legend('Model','Data');
                    legend('boxoff');
                end    
                
            end
            set(gcf, 'Color', 'w');
            [rsq rsq_tot] = an.coefDetermination(data, pred);
            [adj_rsq adj_rsq_tot] = an.adjCoefDetermination(data, pred, sum(p_on)-0.01);
        end
        
        
        
        function [rsq rsq_tot] = coefDetermination(an, data, pred)
            %Calculates the coeficiente of determination            
            
            rsq = zeros(size(data, 1),1);
            rsq_tot = 0;
            for i=1:size(data, 1)
                %Compute the residual values as a vector signed numbers:
                yresid = data(i,:) - pred(i,:);
                %Square the residuals and total them obtain the residual sum of squares:
                SSresid = sum(yresid.^2);
                %Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
                SStotal = (length(data(i,:))-1) * var(data(i,:));
                %SStotal = var(data(i,:));
                %Compute R2 using the formula given in the introduction of this topic:
                rsq(i) = 1 - SSresid/SStotal;
            end   
            %rsq
            %mean(rsq)
            %Total
            if(size(data, 1)>1)
                data=reshape(data,size(data,1)*size(data,2),1);
                pred=reshape(pred,size(pred,1)*size(pred,2),1);
                SSresid = sum((data - pred).^2);
                SStotal = (length(data)-1) * var(data);
                rsq_tot = 1 - SSresid/SStotal;
            end
        end   
        
        function [rsq_adj rsq_adj_tot] = adjCoefDetermination(an, data, pred, n_params)
            %Calculates the adjusted coeficiente of determination            
            
            rsq_adj = zeros(size(data, 1),1);
            
            for i=1:size(data, 1)
                %Compute the residual values as a vector signed numbers:
                yresid = data(i,:) - pred(i,:);
                %Square the residuals and total them obtain the residual sum of squares:
                SSresid = sum(yresid.^2);
                %Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
                SStotal = (length(data(i,:))-1) * var(data(i,:));
                %Compute R2 using the formula given in the introduction of this topic:
                rsq_adj(i) = 1 - SSresid/SStotal * (length(data(i,:))-1)/(length(data(i,:))-n_params-1);
            end   
            
            %Total
            data=reshape(data,size(data,1)*size(data,2),1);
            pred=reshape(pred,size(pred,1)*size(pred,2),1);
            SSresid = sum((data - pred).^2);
            SStotal = (length(data)-1) * var(data);
            rsq_adj_tot = 1 - SSresid/SStotal * (length(data)-1)/(length(data)-n_params-1);
        end 
         
            
    end
    
end

