%%%%
%% Inferring TM model (GUI)
%%%%

function TM_GUI()

close all;
clear all;
clc;

global fig_dist
global fig_MAP;
fig_dist = -1;
fig_MAP = -1;

XLSpath = '';
tags_eTM = {'Name', 'D (s)', 'F (s)', 'U', 'F', 'R^2'};
tags_TMfacil = {'Name', 'D (s)', 'F (s)', 'U', 'R^2'};
tags_TM = {'Name', 'D (s)', 'U', 'R^2'};
colw = {150,40,40,40,40,40};

f = figure('Name','Inferring Tsodyks-Markram Model  v0.1','NumberTitle','off','MenuBar','none','ToolBar','none');
movegui(f,'center');

menu_item_method = uimenu(f, 'Label','Method');
menu_item_method_bayesian = uimenu(menu_item_method, 'Label','Bayesian','Callback',@onBayesianPanel);
%menu_item_method_mse = uimenu(menu_item_method, 'Label','MSE','Callback',@onMSEPanel);

%% Bayesian Panel
bayesian_panel = uipanel('Title','Bayesian Panel','FontSize',12,...
             'Position',[0 0 1 1],'Visible','on');
         
bp_options_panel = uipanel(bayesian_panel, 'Title','Options','FontSize',12,...
             'Position',[0.025 .65 .6 .32],'Visible','on');
         
bpop_age_text = uicontrol(bp_options_panel, 'Style','text','String','Age (days): ',...
          'Position',[2,73,70,25]);
      
bpop_age_edit = uicontrol(bp_options_panel, 'Style','edit','String','-1','BackgroundColor','w',...
          'Position',[65,77,40,25]);
         
%bpop_recpulse_cbox = uicontrol(bp_options_panel, 'Style','checkbox','String','Recovery pulses?',...
%          'Position',[10,55,120,25],...
%          'Callback',{@onRecPulses});
      
bpop_minpulses_text = uicontrol(bp_options_panel, 'Style','text','String','Min pulses: ',...
          'Position',[140,50,70,25]);
      
bpop_minpulses_edit = uicontrol(bp_options_panel, 'Style','edit','String','5','BackgroundColor','w',...
          'Position',[140+65,55,40,25]);%,'Max', 150, 'Min', 0.1);
      
bpop_freq_text = uicontrol(bp_options_panel, 'Style','text','String','Freq (Hz): ',...
          'Position',[2,25,70,25]);
      
bpop_freq_edit = uicontrol(bp_options_panel, 'Style','edit','String','30','BackgroundColor','w',...
          'Position',[65,31,40,25]);%,'Max', 150, 'Min', 0.1);
      
bpop_model_text = uicontrol(bp_options_panel, 'Style','text','String','Model: ',...
          'Position',[140,16,70,25]);
      
bpop_model_popup = uicontrol(bp_options_panel, 'Style','popupmenu','String',{'Extended TM (4 params)', 'TM with facil. (3 params)', 'TM (2 params)'},'BackgroundColor','w',...
          'Position',[140+60,20,120,25],...
          'Callback',{@updateModel});
         
bp_table_panel = uipanel(bayesian_panel, 'Title','MAP Solutions','FontSize',12,...
             'Position',[0.025 .1 .9 .40],'Visible','on');         

bp_table_panel = uitable(bp_table_panel,'ColumnName',tags_eTM,'ColumnWidth',colw,...
             'Position',[20 20 460 100]);
         
%bp_plot = uipanel(bayesian_panel, 'Title','Plot','FontSize',12,'Position',[0 0 1 1],'Visible','on');         


bp_loadfile_pbutton = uicontrol(bayesian_panel, 'Style','pushbutton','String','Load XLS File',...
          'Position',[400,345,80,25],...
          'Callback',@(obj,e)getSTP_XLSFile);
      
bp_loadfile_text = uicontrol(bayesian_panel, 'Style','text','String','',...
          'Position',[370,345-25,140,25]);
      
bp_run_pbutton = uicontrol(bayesian_panel, 'Style','pushbutton','String','Run',...
          'Position',[400,280,80,25],'Enable','off',...
          'Callback',{@runBayesianInference});
      
bp_save_pbutton = uicontrol(bayesian_panel, 'Style','pushbutton','String','Save plots',...
          'Position',[400,240,80,25],'Enable','off',...
          'Callback',{@savePlots});      
      
    
         
%% MSE Panel

%mse_panel = uipanel('Title','MSE Panel','FontSize',12,...
%             'Position',[0 0 1 1],'Visible','off');
         
%msep_loadfile_pbutton = uicontrol(mse_panel, 'Style','pushbutton','String','Load XLS File',...
%          'Position',[.5,.5,70,25],...
%          'Callback',{@getSTP_XLSFile});         
        

function file_path=getSTP_XLSFile(hObject, eventdata)   
    
    [FileName,PathName,FilterIndex] = uigetfile('*.xls','Select the STP XLS file', ['data' filesep() 'ePhys']);
    if(FileName~=0)
        set(bp_loadfile_text,'String',FileName);
        set(bp_run_pbutton,'Enable','on');        
        XLSpath = [PathName FileName];
    end
end

function runBayesianInference(hObject, eventdata)   
    set(bp_run_pbutton,'Enable','Off','String','Running..');
    if(fig_dist>-1 && fig_MAP>-1)
        set(bp_save_pbutton,'Enable','off');
        close(fig_dist);
        close(fig_MAP);
    end
    pause(1);
    [sols names rsq fig_dist fig_MAP] = TM_Bayesian(XLSpath, get(bpop_model_popup,'Value'),...
        str2double(get(bpop_freq_edit,'String')), str2double(get(bpop_minpulses_edit,'String')), str2double(get(bpop_age_edit,'String')));
    
    if(~isempty(sols))
        sols = [names num2cell([sols rsq])];
        set(bp_table_panel,'Data',sols);    
        set(bp_save_pbutton,'Enable','on');
        disp('>>> Done')
    end
    set(bp_run_pbutton,'Enable','On','String','Run');
end


function savePlots(hObject, eventdata)   
    [file,path] = uiputfile('*.pdf','Save Plots As','TM_plot');
    if(isequal(file,0)==0 && isequal(file,0)==0)
        export_fig([path file(1:end-4) '_posterior.pdf'], fig_dist);
        export_fig([path  file(1:end-4) '_MAP.pdf'], fig_MAP);
        disp([path file(1:end-4) '_posterior.pdf']);
        disp([path file(1:end-4) '_MAP.pdf']);
    end    
end

function updateModel(hObject, eventdata)
    v=get(bpop_model_popup,'Value');
    if(v==1)
        set(bp_table_panel,'ColumnName',tags_eTM);
    elseif (v==2)
        set(bp_table_panel,'ColumnName',tags_TMfacil);
    elseif (v==3)
        set(bp_table_panel,'ColumnName',tags_TM);
    end
    
end

function onRecPulses(hObject, eventdata)
    s=get(bpop_recpulse_cbox,'Value');
    if(s)
        set(bpop_minpulses_edit,'Enable','off');
        set(bpop_minpulses_edit,'String','5');
    else
        set(bpop_minpulses_edit,'Enable','On');
    end    
    
end

function onBayesianPanel(hObject, eventdata)   
    set(bayesian_panel,'Visible','on');
    set(mse_panel,'Visible','off');
end

function onMSEPanel(hObject, eventdata)   
    set(bayesian_panel,'Visible','off');
    set(mse_panel,'Visible','on');
end

end