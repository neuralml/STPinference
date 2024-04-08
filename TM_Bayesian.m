%
% Bayesian Inference of Tsodyks-Markram Short-term Synaptic Plasticity model
%
%  data_path: path to the xls data file
%  model_version: specifies the TM model formulation (1 - eTM, 2 - TM with facil, and 3 - TM)
%  f: protocol frequency
%  minp: minimum number of pulses
%  age: age to be selected from the .xls file

% TM_Bayesian("data/examplePCPC.xls", 1, 30, 5, -1)

function [ps_map names rsq fig_dist fig_MAP] = TM_Bayesian(data_path, model_version, f, minp, age)

    global dt;
    dt = 1e-3;
    synbox_path = pwd;
    addpath(genpath(synbox_path));
    
    %% 0. Options
    % 0.1 Data
    freqs = f; %Presynaptic pulse frequency
    min_pulses = minp; %Number of pulses
    remove_drugs = {''};
    remove_areas = {''};
    remove_quality = {'Bad'};
    stime = min_pulses/freqs(1)-dt*2;
    
    cconv = 0;
    cnconv = 0;
    fig_dist = -1;
    fig_MAP = -1;
    dens = {};
    
    if(model_version==1)
        tags = {'Depression Timeconstant, D(s)', 'Facilitation Timeconstant, F(s)', 'Release Probability, U', 'Facilitation Rate, f'};
    elseif(model_version==2)
        tags = {'Depression Timeconstant, D(s)', 'Facilitation Timeconstant, F(s)', 'Release Probability, U'};
    elseif(model_version==3)
        tags = {'Depression Timeconstant, D(s)', 'Release Probability, U'};
    end    

    % 0.2 Optimization
    runs = 3;
    p_limits = [1.00e-04  5.00e-04   2;
               1.00e-04   5.00e-04   2;
               1.00e-04   5.00e-04   1;
               1.00e-04   5.00e-04   1;
               1.00e-04   5.00e-04   1];


    if(model_version==1)
        p_on = [1 1 1 0 1]; %Params: DFUf
    elseif(model_version==2)
        p_on = [1 1 1 0 0]; %Params: DFU
    else
        p_on = [1 0 1 0 0]; %Params: DU    
    end
    
    online_plot = 0;
    cond_reuse = 0;
    conv = 1; %Test convergence?

    %% 1. Load Model
    model = Plasticity.STP.Pheno.MarkramTsodyks98();

    %% 2. Specify protocol and data
    condition = Fitting.Loading.ePhys.STP.DL_STP.NONE_COND_NOTNORM;
    
    dl = Fitting.Loading.ePhys.STP.DL_STP(condition);
    dl.freqs = freqs; dl.min_pulses = min_pulses; dl.remove_drugs = remove_drugs; dl.remove_areas = remove_areas; dl.remove_quality = remove_quality;
    dl.age = age;


    [data CVs names STDs] = dl.load(data_path); %Load data
    
    [model.spikes model.stimes] = dl.setInput(stime); %Set spike train
    model.setDttimes(); %Create dttimes in the model (used for a more optimized)
    
    if(isempty(data))
        warning('The options selected do not match any datapoint');
        ps_map = []; names = []; rsq = [];
        return;
    end

    %dl.stds = bsxfun(@times, ones(size(data)), CVs);
    dl.stds = STDs;

    %% 3. Load optimizer
    opt = Fitting.Optimization.MCMC.SliceSampling();
    opt.burnin = 2500;
    opt.n_samples = 7500;
    if(model_version==1)
        opt.width = [2 2 1 1];
    elseif(model_version==2)
        opt.width = [2 2 1];
    else
        opt.width = [2 1];
    end
    opt.thin = 1;
    opt.eval = @Fitting.Loading.ePhys.logpdf;

    model.run_fun = @model.run4Opt_Fast; %Analytical ver.
    
    opt.runs = runs;
    opt.p_limits = p_limits;
    opt.online_plot = online_plot;
    opt.cond_reuse = cond_reuse; % Starts the opt with the best from the other cond

    %% 4. Run Sampling
    opt.p_on = p_on;
    disp(' ');
    disp('>>> Starting posterior sampling..');
    keyboard
    [ps es] = opt.run(data, dl, model, [], tags, []);
    keyboard
    %Obtain Max a posterior params
    disp(' ');
    disp('>>> Getting maximum a posterior (MAP) solutions..');
    ps_map = opt.map_bestsample(opt.eval, ps, p_on, data, model, p_limits, dl.stds, Inf)';
    
    for i=1:size(data,1)
        dens{i} = kde(ps{i}(:,:)', 'rot');
    end

    %Get R
    if(conv && runs>1) %Convergence diagnostic
        disp(' ');
        disp('>>> Testing convergence..')
        C = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
        for ds=1:size(opt.dist,1) %For each dataset
            PSRF = opt.convergence(ds, opt.dist, C, tags, 0);
            PSRF = mean(PSRF,1);
            
            if(sum(PSRF<1.1 & PSRF>0.9)==sum(p_on))                
                disp([names{ds} ' converged!']);
                cconv = cconv + 1;
            else
                disp([names{ds} ' didnt converged!']);
                cnconv = cnconv + 1;
            end    
            %PSRF
        end
        disp([num2str(cconv) '/' num2str(cnconv+cconv) ' data points converged!']);
    end    

    %% 5. Load analyzer
    ps_map = ps_map';
    an = Fitting.Analyser();
    
    disp(' ');
    disp('>>> Plotting..')
    [results ppr rsq rsq_tot adj_rsq adj_rsq_tot] = an.compareResults(data, names, ps_map, model, p_on, 0);
    fig_MAP = gcf;
    
    an.distPlot(tags, p_limits, dens, ps_map);
    fig_dist = gcf;
end