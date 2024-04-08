function [spikes, stimes] = inhreg(t, f, d)
% function [spikes, stimes isi] = inhreg(t, dt, f)
% Inhomogenous regular distributed ISIs
% t - time
% dt - time step
% f - instantaneous rate vector (per timebase)
% Basic assumptions are:
% (1) constant rate (frequency) over a time step
% (2) only a single arrival possible in a time step
% (so time step should be small relative to the rate of change in
% frequency and arrival rate)

global dt;

spikes=zeros(1,round(t/dt));    
%spikes(1)=0;
%tp=0;   % index of previous spike

% for i=2:length(t)
%     cisi=1/f(i);    % current ISI (seconds)
%     if (i-tp)*dt >= cisi
%         spikes(i)=1;
%         tp=i;       % index of new spike
%     end;
% end;

spikes(1:round(1/f(1)/dt):end)=1;
delay=zeros(1,round(d/dt));  
spikes = [delay spikes];

stimes=find(spikes==1)*dt;
%isi=stimes(2:length(stimes))-stimes(1:length(stimes)-1);
