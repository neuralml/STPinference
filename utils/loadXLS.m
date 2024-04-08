function [num txt]=loadXLS(path, dl)
% path: string with the datafile path
% dl: dataloading class

%     total_pulses = dl.pulses
%     steady_state = dl.pulses;
%     n_pulses = total_pulses;
%     remove_quality = dl.remove_quality; %Poor, Ok, Good
%     remove_areas = dl.remove_areas; %e.g. SC >> Todo if needed
%     remove_drugs = dl.remove_drugs; %e.g. Ro
%     min_pulses = 0; %Use min number of pulses
%     age_freqs = [10 30]; %Frequencies to use in the age-dependent plot

    [num txt] = xlsread(path);
    %[num txt raw] = xlsread('/Users/ruicosta/Documents/UoE/PhD/own/experiments/experimental/analysis/data_before_after_freq.xls');
    
    %ages = num(:,16);
    
    
    
%     %Find how many diff frequencies
%     num(num(:,1)==30, 1) = 15;
%     freqs = sort(unique(num(:,1)));
%     n_freq = length(freqs);
%     
% 
%     %Create array to store the means and std
%     results_before = zeros(2,n_freq,n_pulses);
%     results_after = zeros(2,n_freq,n_pulses);
%     results_ppr_b = zeros(2,n_freq,n_pulses-1);
%     results_ppr_a = zeros(2,n_freq,n_pulses-1);
%     
%     num_before = [];
%     num_after = [];
%     ppr_b = [];
%     ppr_a = [];
%     
%     %Normalize
%     c1 = 1;
%     c2 = 1;
%     norm = 0;
%     
%     debug = {};
%     for a=1:size(num,1)
%         %num(a,2:1+total_pulses)
%         if(mod(a,2)==1)%If even, then 1,3,5,7 (before wash-in)
%             norm = num(a,2); %1st EPSP (to normalize)
%             num_before(c1,:) = [num(a,1) num(a,2:1+total_pulses)./norm];
%             for p=3:1+n_pulses%Calculates the PPR the 5 pulses
%                 ppr_b(c1,p-2) = (num(a,p)-num(a,2))/num(a,2);
%             end
%             c1=c1+1;
%         else% (after wash-in)
%             num_after(c2,:) = [num(a,1) num(a,2:1+total_pulses)./norm];
%             if(num(a,1)==30)
%                 disp([num2str(num(a,2)) '/' num2str(norm)]);
%             end
%             %num_after(c2,1+steady_state) = [num(i,1+steady_state)./num(i,2)]; %normalize last pulse by the 1st EPSP after
%             
%             debug{c2} = [num2str(num(a,1)) [txt(a,1) txt(a,5)] mat2str(num_after(c2,1:1+total_pulses))];
%             for p=3:1+n_pulses%Calculates the PPR for the 5 pulses
%                 ppr_a(c2,p-2) = (num(a,p)-num(a,2))/num(a,2);
%             end
%             c2=c2+1;            
%         end    
%         
%     end
%     
%     debug{:,:}
%    
%     
% %     num_before(1:10,1:6)
% %     num_after(1:10,1:6)
% %     pause
%     
%     %Get data Before washin
%     for i=1:n_freq
% %         freqs(i)
%         numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%         
%         results_before(1,i,:) = calc_mean(num_before(numi, 2:1+n_pulses)); %Mean
%         results_before(2,i,:) = calc_std(num_before(numi, 2:1+n_pulses))./sqrt(sum(numi)); %SEM
%         results_ppr_b(1,i,:) = calc_mean(ppr_b(numi,:)); %Mean
%         results_ppr_b(2,i,:) = calc_std(ppr_b(numi,:))./sqrt(sum(numi)); %SEM
%     end    
%     
%     
%     %Get data After washin
%     for i=1:n_freq
% %         freqs(i)
%         numi = ((num_after(:,1) == freqs(i)) & ~isnan(num_after(:,2)));
%         
%         results_after(1,i,:) = calc_mean(num_after(numi, 2:1+n_pulses)); %Mean
%         results_after(2,i,:) = calc_std(num_after(numi, 2:1+n_pulses))./sqrt(sum(numi)); %SEM
%         results_ppr_a(1,i,:) = calc_mean(ppr_a(numi,:)); %Mean
%         results_ppr_a(2,i,:) = calc_std(ppr_a(numi,:))./sqrt(sum(numi)); %SEM
%     end    
%     
%     g = [0.6 0.6 0.6];
%     
%     freqs_tmp=freqs;
%     %% Plot
%     if(length(freqs)==1)
%         figure;
%         r=reshape(results_after(1,:,1), n_freq, 1)*100;
% 
%         bar(freqs, r, 1, 'b', 'EdgeColor', [1 1 1], 'BarWidth', 0.5); %1st EPSP (Mean)
%         xlim([freqs(1)-1 freqs(1)+1]);
%         ylim([0 110]);
%         hold on;
%         errorbar(freqs,results_after(1,:,1)*100,results_after(2,:,1)*100,'-o','Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1], 'LineWidth', 1.1); %1st EPSP (STD)
%         hold on;
%         %scatter plot
%         for i=1:n_freq
%             numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%             text(freqs(i)-0.025, 10, ['n=' num2str(sum(numi))],'Color',[1 1 1]);
%             scatter(freqs(i).*ones(1,sum(numi)),num_after(numi, 2)*100,[],[0.5 0.5 1],'MarkerFaceColor','w');
%         end    
%         ylabel('EPSP_{After}^{1}/ EPSP_{Before}^{1} (%)');
%         legend('After Washin');
%         %set(gca,'XTick', [freqs(freqs<7.5); 7.5; freqs(freqs>=7.5)]);
%         set(gca,'XTick', freqs);
%         set(gca, 'Box', 'off' );
%         %xlabel('Frequency (Hz)');
%         %set(gca,'XTickLable', 'After');
%         
%     else
%         
%         figure;
%         subplot(3,1,1);
%         %results_after(1,:,1)
%         %r=reshape(results_after(1,:,1:5), n_freq, 5)*100
%         r=reshape(results_after(1,:,1), n_freq, 1)*100;
% 
%         %bar3(freqs, r, 0.5); %1st EPSP (Mean)
%         %bar(freqs, r, 3.5, 'grouped', 'EdgeColor', [1 1 1]); %1st EPSP (Mean)
%         bar(freqs, r, 1, 'b', 'EdgeColor', [1 1 1]); %1st EPSP (Mean)
%         hold on;
%         if(nargin>1)
%             plot(0,0,'--g', 'Color', [0 0.8 0], 'LineWidth', 1.1);
%             hold on;
%         end    
%         xlim([-0.7 16]);
%         hold on;
%         errorbar(freqs,results_after(1,:,1)*100,results_after(2,:,1)*100,'-o','Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1], 'LineWidth', 1.1); %1st EPSP (STD)
%         hold on;
%         %scatter plot
%         for i=1:n_freq
%             numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%             text(freqs(i)-0.33, 10, ['n=' num2str(sum(numi))],'Color',[1 1 1]);
%             scatter(freqs(i).*ones(1,sum(numi)),num_after(numi, 2)*100,[],[0.5 0.5 1],'MarkerFaceColor','w');
%         end    
%         ylabel('EPSP_{After}^{1}/ EPSP_{Before}^{1} (%)');
%         
%         if(nargin<2)
%             legend('After Washin');
%         else
%             hold on;
%             plot(model_freqs, model_results, '--g', 'Color', [0 0.8 0], 'LineWidth', 1.1);
%             legend('After Washin', 'Model Prediction');
%         end
%         
%         
%         %set(gca,'XTick', [freqs(freqs<7.5); 7.5; freqs(freqs>=7.5)]);
%         set(gca,'XTick', freqs);
%         set(gca, 'Box', 'off' );
%         freqs_tmp(end) = 30;
%         set(gca,'XTickLabel',freqs_tmp);
%         
%         hcolors = hot(10); %Use hot colormap
% 
%         fitSigmoid();
%         
% 
%         %% 1st PPR
%         freqs(1) = []; %Remove 0.1Hz
%         results_ppr_b(:,1,:) = [];
%         results_ppr_a(:,1,:) = [];
%         n_freq = n_freq-1;
%         
%         h=subplot(3,1,2);
%         r=reshape([results_ppr_b(1,:,1) results_ppr_a(1,:,1)], n_freq, 2)*100;
%         map = [1 0 0; 0 0 1];
%         colormap(h,map);
%         
%         bar(freqs, r, 2.2, 'group', 'EdgeColor', [1 1 1]); %1st EPSP (Mean)
%         hold on;
%         plot(freqs, r(:,2)-r(:,1), '-o', 'Color', hcolors(1,:), 'LineWidth', 1.1,'MarkerFaceColor','w');
%         hold on;
%         xlim([-0.7 16]);
%         errorbar(freqs-0.14,results_ppr_b(1,:,1)*100, results_ppr_b(2,:,1)*100,'-o','Color',[1 0.5 0.5],'MarkerFaceColor',[1 0.5 0.5], 'LineWidth', 1.1); %1st EPSP (STD)
%         hold on;
%         errorbar(freqs+0.12,results_ppr_a(1,:,1)*100, results_ppr_a(2,:,1)*100,'-o','Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1], 'LineWidth', 1.1); %1st EPSP (STD)
%         hold on;
%         for i=1:n_freq
%             numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%             %text(freqs(i)-0.3, 10, ['n=' num2str(sum(numi))],'Color',[1 1 1]);
%             scatter(freqs(i).*ones(1,sum(numi))+0.12,ppr_a(numi)*100,[],[0.5 0.5 1],'MarkerFaceColor','w');
%             scatter(freqs(i).*ones(1,sum(numi))-0.14,ppr_b(numi)*100,[],[1 0.5 0.5],'MarkerFaceColor','w');
%         end
%         legend('Before Washin','After Washin','PPR_{After} - PPR_{Before}');
%         
%         ylabel('PPR_1 (%)');
%         set(gca,'XTick',freqs);
%         set(gca, 'Box', 'off' );
%         set(gca,'XTickLabel',freqs_tmp);
% 
%         %(EPSP_n-EPSP_1)/EPSP_1
%         h=subplot(3,1,3);
%         %bar3(freqs, r, 0.5); %1st EPSP (Mean)
%         %bar(freqs, r, 3.5, 'grouped', 'EdgeColor', [1 1 1]); %1st EPSP (Mean)
%         %map = [1 0 0; 0 0 1];
% 
%         %bar(freqs, r, 2.2, 'group', 'EdgeColor', [1 1 1]); %1st EPSP (Mean)
%         xlim([-0.7 31]);
%         %hold on;
%         l = {};
%         for i=1:n_pulses-1
%             r=reshape([results_ppr_b(1,:,i) results_ppr_a(1,:,i)], n_freq, 2)*100;
%             plot(freqs, r(:,2)-r(:,1), '-o', 'LineWidth', 1.1, 'Color', hcolors(i,:),'MarkerFaceColor','w');
%             l{i} = ['PPR_' num2str(i)];
%             hold on;
%         end    
%         %hold off;
% 
%     %     hold on;
%     %     errorbar(freqs-0.14,results_before(1,:,5)*100, results_before(2,:,5)*100,'-o','Color',[1 0.5 0.5],'MarkerFaceColor',[1 0.5 0.5]); %1st EPSP (STD)
%     %     hold on;
%     %     errorbar(freqs+0.12,results_after(1,:,5)*100, results_after(2,:,5)*100,'-o','Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1]); %1st EPSP (STD)  
%     %     hold on;
%     %     for i=1:n_freq
%     %         numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%     %         %text(freqs(i)-0.3, 10, ['n=' num2str(sum(numi))],'Color',[1 1 1]);
%     %         scatter(freqs(i).*ones(1,sum(numi))+0.12,num_after(numi, 6)*100,[],[0.5 0.5 1],'MarkerFaceColor','w');
%     %         scatter(freqs(i).*ones(1,sum(numi))-0.14,num_before(numi, 6)*100,[],[1 0.5 0.5],'MarkerFaceColor','w');
%     %     end
%         %legend('Before Washin','After Washin','Increase');
%         legend(l);
%         ylabel('PPR_{After} - PPR_{Before} (%)');
% 
%         set(gca,'XTick', freqs);
%         freqs_tmp
%         pause
%         set(gca,'XTickLabel',freqs_tmp);
% 
%         set(gca, 'Box', 'off' );
%         xlabel('Frequency (Hz)');
% 
%         %% Plot Age-dependence
%         %plotAge()
%         %set(gca, 'Box', 'off' );
%         %num
%         %txt
%     end
%     
%     
%     
%     function mean_vector = calc_mean(m)
%         %Compute mean by columns
%         mean_vector = zeros(1,size(m,2));
%         
%         for x=1:size(m,2)
%             c = 0;
%             t = 0;
%             for j=1:size(m,1)
%                 if ~isnan(m(j,x))
%                     c = c+1;
%                     t = t+m(j,x); 
%                 end    
%             end
%             mean_vector(x) = t/c;
%         end
%         %Replace NaN by copies of the last
%         aux=find(~isnan(mean_vector));
%         mean_vector(isnan(mean_vector)) = mean_vector(aux(end));
%         
%     end
% 
%     function std_vector = calc_std(m)
%         %Compute std by columns
%         std_vector = zeros(1,size(m,2));
%         for x=1:size(m,2)
%             t = [];
%             c=0;
%             for j=1:size(m,1)
%                 if ~isnan(m(j,x))
%                     c = c+1;
%                     t(c) = m(j,x); 
%                 end    
%             end
%             std_vector(x) = std(t);
%         end
%         
%         aux=find(~isnan(std_vector));
%         std_vector(isnan(std_vector)) = std_vector(aux(end));
%         
%     end
% 
%     function fitSigmoid()
%         %% Fit a sigmoid to data
%         %%
%         xf = [];
%         %y = 1./(5 + 10./(1+exp(-(x-40)/10)) + randn(size(x)));
%         yf = [];
% 
%         %Build y
%         for i=1:n_freq
%             numi = ((num_before(:,1) == freqs(i)) & ~isnan(num_before(:,2)));
%             yf = [yf; (num_after(numi, 2)*100)];        
%             xf = [xf; (freqs(i).*ones(1,sum(numi)))']; 
% 
%             %scatter(freqs(i).*ones(1,sum(numi)), num_after(numi, 2)*100,[],[0.5 0.5 1],'MarkerFaceColor','w');
%         end 
%         %xf
%         %yf
% 
%         %Add hypothetical data
% %         yf = [yf; 65; 85; 70];
% %         xf = [xf; 7.5; 6.5; 8.5];
%         %yf = [yf; 75; 97];
%         %xf = [xf; 10; 0.1];
% 
%         %figure;
%         %plot(xf,yf,'bo');
%         %hold on;
% 
%         %%
%          fun = @(p,x) (p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4))));
%     %     f = @(p,x) (p(1) ./ (1 + exp(-(x-p(2)))));
%          %[p,resid,J,Sigma] = nlinfit(xf,yf,fun,[67.78 24.32 7.771 -0.7154]);
%          %p
%          %line(xf, fun(p,xf), 'color', 'k')
%     %      ci = nlparci(p,resid,'jacobian',J)
%     %     line(xf, f(ci(:,1),xf), 'color', 'g')
%     %     line(xf, f(ci(:,2),xf), 'color', 'b')
% 
% 
%         s = fitoptions('Method','NonlinearLeastSquares',...
%                    'Display','iter', 'StartPoint', [65.73 33.14 8.37 -0.8583]);
% 
%         f = fittype('a + b ./ (1 + exp(-(x-m)/s))', 'options', s);       
% 
%         [ft gof output] = fit(xf, yf, f); 
%         ft
%         gof
%         category(ft)
%         
%         %fdata = feval(ft,xf); 
%         %I = abs(fdata - yf) > 1.5*std(yf); 
%         %outliers = excludedata(xf,yf,'indices',I);
% 
%         %[ft2 gof output] = fit(xf, yf, f, 'Exclude', outliers);
%         
%         %ci = nlparci([ft.a ft.b ft.m ft.s], output.residuals, 'jacobian',output.Jacobian)
%         %[ft gof]=fit(xf,yf,'a + b ./ (1 + exp(-(x-m)/s))','start',[0 20 50 5]) 
% 
%         res = fun([ft.a ft.b ft.m ft.s], 0.1:0.01:30);
%         plot(0.1:0.01:30, res, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.05);
%         %plot(ft)
%         hold on;
%         
%         confs = confint(ft,0.95)
%         gof_d = 0.1; % 5% change
%         %[ft.a-ft.a*gof_d ft.b-ft.b*gof_d ft.m-ft.m*gof_d ft.s-ft.s*gof_d]
%         less = fun([ft.a-ft.a*gof_d ft.b-ft.b*gof_d ft.m-ft.m*gof_d ft.s-ft.s*gof_d], 0.1:0.01:30)';
%         more = fun([ft.a+ft.a*gof_d ft.b+ft.b*gof_d ft.m+ft.m*gof_d ft.s+ft.s*gof_d], 0.1:0.01:30)';
%         %confs(2,4) = ft.s+ft.s*gof_d;
%         %less = fun([confs(1,:)], 0.1:0.01:30)';
%         %more = fun([confs(2,:)], 0.1:0.01:30)';
%         
% %         plot(0.1:0.01:30, less, '-g', 'LineWidth', 1.05);
% %         hold on
% %         plot(0.1:0.01:30, more, '-k', 'LineWidth', 1.05);
% %         hold on
% 
%         
% %         less2 = fun([confs(1,1) confs(1,2) ft.m+confs(1,3) ft.s+confs(1,4)], 0.1:0.01:30)';
% %         more2 = fun([confs(2,1) confs(2,2) ft.m+confs(2,3) ft.s+confs(2,4)], 0.1:0.01:30)';
%         %plot(xf, bounds, 'm--');
%         
%         %xf
%         %[less more-less]
%         h = area(0.1:0.01:30, [less more-less]);
% 
%         %set(gca,'Layer','top');    
%         set(h(1),'FaceColor', [1 1 1])
%         set(h(2),'FaceColor', [0.8 0.8 0.8])
%         set(h,'EdgeColor', [1 1 1]);
% 
%         handles = get(gca,'children');
%         handles = [handles; handles(3); handles(1:2)];
%         handles(1:3) = [];
%         set(gca,'children', handles);
% %         figure
% %         plot(ft,'r-',xf,yf,'k.',outliers,'m*') 
% %         hold on
% %         plot(ft2,'c--') 
% % 
% %         figure
%     end    
% 
% 
%     function plotAge()
%         figure;
%         
%         sp=subplot(1,2,1);
%         %r=reshape(results_after(1,:,1), n_freq, 1)*100;
%         ages(1:2:end) = [];
%         ages
%         c=0;
%         xx = [];
%         yy = [];
%         for j=1:length(age_freqs)
%             numi = ((num_before(:,1) == age_freqs(j)) & ~isnan(num_before(:,2)))
%             %text(freqs(i)-0.33, 10, ['n=' num2str(sum(numi))],'Color',[1 1 1]);
%             
%             xx = [xx; ages(numi)];
%             yy = [yy; num_after(numi, 2)*100];
%             c=c+sum(numi);
%         end
%         scatter(xx, yy,[],[0.5 0.5 1],'MarkerFaceColor','k');
%         hold on;
%         c
%         
%         % Find x values for plotting the fit based on xlim
%         axesLimits1 = xlim(sp);
%         xplot1 = linspace(axesLimits1(1), axesLimits1(2));
% 
%         % Find coefficients for polynomial (order = 1)
%         fitResults1 = polyfit(xx, yy, 1)
%         % Evaluate polynomial
%         yplot1 = polyval(fitResults1, xplot1)
%         % Plot the fit
%         [cc pp] = corrcoef(xx,yy)
%         plot(xplot1, yplot1, '-r', 'DisplayName',['r=' num2str(cc(1,2)) ' p=' num2str(pp(1,2))]);
%         
%         ylabel('EPSP_{After}^{1}/ EPSP_{Before}^{1} (%)');
%         xlabel('Age (days)');
%         legend(sp,'show');
%         
%         %set(gca,'XTick', freqs);
%     
%         sp=subplot(1,2,2);
%         xx = [];
%         yy = [];
%         for j=1:length(age_freqs)
%             numi = ((num_before(:,1) == age_freqs(j)) & ~isnan(num_before(:,2)));
%             xx = [xx; ages(numi)];
%             yy = [yy; ppr_a(numi)*100-ppr_b(numi)*100];
%         end
%         scatter(xx, yy,[],[0.5 0.5 1], 'MarkerFaceColor','k');
%         hold on;
%         
%         % Find x values for plotting the fit based on xlim
%         axesLimits1 = xlim(sp);
%         xplot1 = linspace(axesLimits1(1), axesLimits1(2));
% 
%         % Find coefficients for polynomial (order = 1)
%         fitResults1 = polyfit(xx, yy, 1)
%         % Evaluate polynomial
%         yplot1 = polyval(fitResults1, xplot1)
%         % Plot the fit
%         [cc pp] = corrcoef(xx,yy)
%         plot(xplot1, yplot1, '-r', 'DisplayName',['r=' num2str(cc(1,2)) ' p=' num2str(pp(1,2))]);
%         
%         
%         ylabel('PPR_{After} - PPR_{Before}');
%         xlabel('Age (days)');
%         legend(sp,'show');
%     end    
end