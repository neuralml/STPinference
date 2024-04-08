classdef MarkramTsodyks98 < Plasticity.LearningRule
    %STP model by Markram et al. PNAS (1998) vol. 95 (9) pp. 5323-8
    % This model considers facilitation and depression
    
    
    properties
         tau_r = 130e-3;
         tau_p = 530e-3;
         tau_in = 10e-3;
         
         freq = 30; %Hz
         
         p_b = 0.03;
         tau_mem = 60e-3;
         n_b = 1;
         f = 0.06;
         v_rest = -65e-3;
         R_in = 1;
         Ase = 1;
         E = 0;
         s;
         n;
         p;
         v;
         
         
         %Storage variables
         nA;
         pA;
         vA;
         eA;
         IA;
         I_peaks;
         verbose = 0;
         
         silentAS_flag = 0;
         
         noise=0;
         
         Isyn = 0;
         
         dt;
         spikes; %Stores spike train
         stimes; %Stores the isis
         dttimes; %Stores the diff between the spike times
         
         syn_input;         
    end
    
    
    methods
        %
        %Contructor
        %
        function stp = MarkramTsodyks98(tau_r, tau_p, tau_in, Pb, tau_mem, facil, Ase)
            
            disp('-- Created MarkramTsodyks98 STP Model --');
            
            if  nargin > 0
                stp.tau_r = tau_r;
                stp.tau_p = tau_p;
                stp.tau_in = tau_in;
                stp.tau_mem = tau_mem;
                
                stp.Ase = Ase;
                stp.p_b = Pb;
                stp.f = facil;
            end    
            stp.n = stp.n_b;
            stp.v = stp.v_rest;
            stp.p = stp.p_b;
            
            global dt;
            stp.dt = dt;
        end
        
        
        function [peaks] = run4Opt_Fast(stp)
            %Run the learning rule for stime milliseconds with one spike train
            %Fast version: Uses the DE solution (exponential decay)
            %in: spike train (0,1)
            
            
            peaks = zeros(1,length(stp.stimes));
            
            p_old = stp.p_b;
            n_old = stp.n_b;
            
            for t=1:length(stp.dttimes) %Simulation cycle
                
                n_new = stp.n_b - (1 - n_old*(1-p_old))*exp(-stp.dttimes(t)/stp.tau_r);
                p_new = stp.p_b + exp(-stp.dttimes(t)/stp.tau_p)*(p_old + stp.f*(1-p_old) - stp.p_b);
                                
                
                %[p_new n_new]
                peaks(t) = stp.Ase*p_old*n_old;
                p_old = p_new;
                n_old = n_new;
            end
            peaks(t+1) = stp.Ase*p_old*n_old;
        end
        
        function reset(stp)
            stp.n = stp.n_b;
            stp.v = stp.v_rest;
            stp.p = stp.p_b;
            stp.E = 0;
        end
        
        function setParams(stp, p, p_on)
            x = [stp.tau_r stp.tau_p stp.p_b stp.Ase stp.f];
            x(p_on==1) = p;
            
            xCell = num2cell(x);
            [stp.tau_r stp.tau_p stp.p_b stp.Ase stp.f] = xCell{:};            
            
            
            
            if(length(p_on)>4) % (3 params)
                if(p_on(5)==0) %Set stp.f = stp.p_b if stp.f is not being optimized
                   stp.f = stp.p_b;
                end    
            end    
            
            if(p_on(2)==0) %Turns off facilitation (2 params)
                stp.f = 0;
                %stp.p_b = 0;
            end
            
        end
        
        function setDttimes(stp)
            stp.dttimes = (stp.stimes(2:end)-stp.stimes(1:end-1));
        end
        
    end
    
end

