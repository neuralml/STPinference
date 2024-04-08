%% Applies MCMC variant slice sampling (uses External procedures)
%

classdef SliceSampling < Fitting.Optimization.MCMC.MCMC
    
    properties
        cond = 0;
        width = 0.5;
        thin = 1;
        burnin = 2000;
        n_samples = 2000;
        dist = []
    end
    
    methods
        
        %
        %Contructor
        %
        function opt = SliceSampling()
            disp('-- Created Slice Sampling --');
            
            opt = opt@Fitting.Optimization.MCMC.MCMC(); % call superclass constructor
        end
        
        
        function [x_bests llh_best] = run(opt, data, dl, models, ps, names, init_cond_opt)
          %e_bests = ones(size(data,1),1).*Inf; %Default best
          llh_best = zeros(size(data,1), 1);
          x_bests = cell(size(data,1), 1);
          opt.dist = zeros(size(data,1), opt.runs, opt.n_samples, sum(opt.p_on==1)); %Stores the samples
          h = 0;
          peaks = 0;
          init_pop = {};
          %an = Fitting.Analyser();
          
          
          if(opt.norm2before)
              peaks = init_cond_opt;
              init_cond_opt = []; 
              ps_peaks = zeros(size(data));
              %Get output for ps_before
              for i=1:size(ps,2)
                models(i).setParams(ps(:,i), opt.p_on);
                models(i).reset();
                [ps_peaks(i,:)] = models(i).run_fun(); %Runs model
              end    
          else
              peaks = ones(size(data,1),1).*Inf;
              ps_peaks = zeros(size(data,1),1);
          end
          
          tic;
          %Run N times
          for i=1:size(data,1)
              %all_params = cell(1,length());
              
              disp(['Sampling - MCMC - Slice Sampling: ' num2str(i) ' out of ' num2str(size(data,1))]);
              
              for l=1:opt.runs
                  disp(['  Run: ' num2str(l)]);
                  
                  if(opt.online_plot)
                      figure;
                      %Plot data
                      scatter(1:length(data(i,:)),data(i,:),'g','filled');
                      hold on;
                      h=scatter(1:length(data(i,:)),data(i,:),'r','filled');
                      hold on;
                  end
                  
                  
                      
                  
                  if(isempty(init_cond_opt))
                    init_cond = rand(1, sum(opt.p_on)).*opt.p_limits(opt.p_on==1,end)'; %Random from Uniform
                  else
                    init_cond = init_cond_opt(i,:);
                  end  
                   
                  if(opt.cond_reuse)%Set the initial condition eq to the previous best                      
                      if(opt.cond)
                        ps_bef = ps(i,:);
                        tmp = opt.p_on;
                        tmp(opt.p_on==0)=[]; %Removes non-used params
                        init_cond = ps_bef(tmp==1);
                      end
                  end 
                  
                  while(1)
                     %try   
                            %keyboard
                          [xt_best] = slicesample(init_cond, opt.n_samples, 'logpdf', @(x) opt.eval(x, data(i,:), models(i), opt.p_on, h, [peaks(i,:); ps_peaks(i,:)], opt.p_limits, dl.stds(i,:)),'width',opt.width,'thin',opt.thin,'burnin',opt.burnin);
                           
                      break;
                     %catch exception
                         %exception
                         %keyboard
                         %fprintf(1,'There was an error! The message was:\n%s', exception.message);
                     %end
                  end
                  
                  opt.dist(i, l, :, :) = xt_best;
                  
                if(opt.online_plot)
                  figure
                  hist(xt_best, 100);
                  legend(names);
                end  


                if(isempty(x_bests{i}))
                    x_bests{i} = xt_best;
                else
                    x_bests{i} = [x_bests{i}; xt_best];
                end
              end
              
              %Trace plot
              if(opt.online_plot)
                opt.plotChains(i, opt.dist, names, models(i), data(i,:), dl.stds(i,:));
              end
          end          
          toc;
          
        end
        
        
        
        function plotChains(opt, data_i, dist, names, model, data, stds)
            figure;
            
            %Calculate the LLH
            llh = zeros(size(dist,2),size(dist,3));
           
            for j=1:size(dist,2)
                for i=1:size(dist,3)
                    aux = reshape(dist(data_i, j, i, :), 1, size(dist,4));
                    llh(j,i) = Fitting.Loading.ePhys.chiSquare(aux, data, model, opt.p_on, [], Inf, opt.p_limits, stds);
                end                
            end    
            
            for i=1:sum(opt.p_on==1)
                subplot(round(sum(opt.p_on==1)/min(sum(opt.p_on==1),2))+1, min(sum(opt.p_on==1),2), i);
                
                aux = reshape(dist(data_i, :, :, i), size(dist,2), size(dist,3));
                plot(1:size(aux,2), aux);
                xlabel('Iterations');
                ylabel('Theta');
                title(names{i});
                ylim([0 opt.p_limits(i,3)]);
            end    
            subplot(round(sum(opt.p_on==1)/min(sum(opt.p_on==1),2))+1, min(sum(opt.p_on==1),2), i+1);  
            plot(1:size(aux,2), min(llh,1));
            ylim([0 1+0.2]);
            xlabel('Iterations');
            ylabel('LLH');
        end
        
    end
    
end

