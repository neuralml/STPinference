%%
%% Generic class for MCMC
%%

classdef MCMC < Fitting.Optimization.Optimization
    
    properties
       nbins = 100;
       conv = 0;
    end
    
    methods
        
        %
        %Contructor
        %
        function opt = MCMC()
            disp('-- Created MCMC --');
            
            opt = opt@Fitting.Optimization.Optimization(); % call superclass constructor
        end
        
        
        function params = run(opt, data)
            
        end 
        
        function [p_final p_L] = map(opt, dens, p_on)
            p_final = zeros(length(dens), sum(p_on==1));
            p_L = zeros(length(dens), sum(p_on==1));

            for i=1:length(dens)                   
                %for j=1:sum(p_on==1) %Get MAP over each marginal
                    p_final(i,:) = max(dens{i})';
                %end
            end
            
        end
        
        function [p_final p_L] = map_normin(opt, dist, p_on)
            p_final = zeros(length(dist), sum(p_on==1));
            p_L = zeros(length(dist), sum(p_on==1));
            npoints = 1000;
            
            for i=1:length(dist)
                for j=1:sum(p_on==1)                    
                    [f,x] = ksdensity(dist{i}(:,j), 'npoints', npoints);
                    aux = x(f==max(max(f)));
                    p_final(i,j) = aux(1)-x(f==min(min(f)));
                    p_L(i,j) = max(max(f));
                    
                end    
            end
        end   
        
        function [p_final L] = map_bestsample(opt, llh, dist, p_on, data, models, p_limits, stds, pre_peaks, ps)
            p_final = zeros(length(dist), sum(p_on==1));
            L = -ones(length(dist), 1).*Inf;
            
            for i=1:length(dist)
                
                if(pre_peaks~=Inf)
                    models(i).setParams(ps(:,i), opt.p_on);
                    models(i).reset();
                    [ps_peaks] = models.run_fun(); %Runs model
                end    
                for j=1:size(dist{i},1) %Over params                
                    if(pre_peaks==Inf)
                        aux = llh(dist{i}(j,:), data(i,:), models(i), p_on, [], Inf, p_limits, stds(i,:));
                    else
                        aux = llh(dist{i}(j,:), data(i,:), models(i), p_on, [], [pre_peaks(i,:); ps_peaks], p_limits, stds(i,:));
                    end
                    if(aux>L(i))
                        p_final(i,:) = dist{i}(j,:);
                        L(i) = aux;
                    end
                end
            end
        end
        
        
        function [p_final L] = map_marginal(opt, llh, dens, p_on, data, model, p_limits, stds, pre_peaks)
            p_final = zeros(length(dens), sum(p_on==1));
            L = ones(length(dens), 1);
            
            for i=1:length(dens)                   
                for j=1:sum(p_on==1) %Get MAP over each marginal
                    p_final(i,j) = max(marginal(dens{i}, j));
                end
                L(i) = llh(p_final(i,:), data(i,:), model, p_on, [], Inf, p_limits, stds(i,:));
            end
        end
        
        
        function [p_final] = map_mean(opt, dist, p_on)
            p_final = zeros(length(dist), sum(p_on==1));
            
            for i=1:length(dist)        
                p_final(i,:) = mean(dist{i}(:,:),1);
            end
        end
        
        function [p_final] = map_median(opt, dist, p_on)
            p_final = zeros(length(dist), sum(p_on==1));
            
            for i=1:length(dist)        
                p_final(i,:) = median(dist{i}(:,:),1);
            end
        end
        
        
        function psrf = convergence(opt, data_i, dist, C, tags, plotOn)
            %% Uses the Gelman and Rubin method
            
            s=size(dist);
            %opt.dist = zeros(size(data,1), opt.runs, opt.n_samples, sum(opt.p_on==1)); %Stores the samples
            m = s(2);            
            n = s(3);
            p = s(4);
            ite=100;
            psrf = zeros(length(ite:ite:n),p); %Potential scale reduction factor
            
            for i=1:length(data_i)
                ni_c = 1;
                for ni = ite:ite:n %For each sampling subset
                    ite2 = ni/2;

                    for pi = 1:p %For each parameter
                        phi = mean(reshape(dist(data_i(i), :, (ni-ite2+1):ni, pi), 1, m*ite2));
                        W_tmp = 0;
                        B_tmp = 0;                    

                        for mi = 1:m %For each chain
                            phi_j = mean(dist(data_i(i), mi, (ni-ite2+1):ni, pi));

                            B_tmp = B_tmp + (phi_j - phi).^2;  %Between chain variance
                            W_tmp = W_tmp + sum((dist(data_i(i), mi, :, pi) - phi_j).^2); %Within chain variance
                        end 

                        %Calculate B and W
                        W = W_tmp/(m*(ite2-1));
                        B = ite2/(m-1)*B_tmp;

                        %Calculate R
                        sigma_est = ((ite2-1)/ite2)*W + B/ite2;
                        %R(pi) = sigma_est/W;
                        psrf(ni_c,pi) = ((m+1)/m)*(sigma_est/W) - (ite2-1)/(m*ite2);
                        %%
                        %disp(['p: ' num2str(pi) ' W: ' num2str(W) ' B: ' num2str(B) ' sigma_est: ' num2str(sigma_est) ' R:' num2str(psrf(pi))]);
                        %%
                    end
                    ni_c = ni_c + 1;
                end
                
                if(plotOn)
                    figure;
                    
                    subplot(round(length(data_i)/min(length(data_i),2))+1, min(length(data_i), 2), i);
                    for pi = 1:p                    
                        plot(ite:ite:n, psrf(:,pi), 'Color' , C(pi,:));
                        hold on;
                    end
                    hold off;
                    ylabel('PSRF');
                    xlabel('Runs');
                    legend(tags);
                end    
            end
        end    
            
    end
    
end

