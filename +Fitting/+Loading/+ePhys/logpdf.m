function [logp peaks] = logpdf(x, data, model, p_on, h, pre_peaks, p_limits, stds)
    %Pre_peaks: data before (1st row) and params before (2nd row)
    
    %keyboard
    
    if(sum(single(x)<single(p_limits(p_on==1,1))') || sum(single(x)>single(p_limits(p_on==1,3))')) %Prior
        logp = -1e+50;
        peaks = 0;
        return;
    end

    %try
    model.setParams(x, p_on);
    model.reset();
    [peaks] = model.run_fun(); %Runs model
    %keyboard
       
    if(pre_peaks(1,1)==Inf) %set data to pre peaks
        peaks = peaks.*sum(data.*peaks./(stds.^2))/sum((peaks.^2)./(stds.^2)); %Scale by A_MAP
    elseif(pre_peaks(1,1)==-1)
        peaks = peaks./peaks(1);
    else        
        peaks = peaks.*sum(pre_peaks(1,:).*(pre_peaks(2,:))./(stds.^2))/sum((pre_peaks(2,:).^2)./(stds.^2)); %Scale by A_MAP
    end

    logp = -sum(1/2.*((peaks(1:end)-data(1:end))./stds).^2);
  end