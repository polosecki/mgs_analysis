close all; clear all

%LIP
%April 11 2013
%LIP
%April 11 2013

% cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
% load(fullfile(cell_file_dir,'LIP_data.mat'));
% {LIP_cell.dir};
%Cells:
%1 --attn effect, weird saccade tuning
%2 -- attn effect--- maybe nothing
%3 -- attn, strange saccade
%4 -- weak attn effect, saccade tuning--weird
%5 -- not clear effect
%6--some attention
%7 -- attention effect
%8 --some attn effect
%9 strong sacc, look at attn, are labels correct? interaction??
%10--
%11 really good, good map, tuning for 135 and 90, slightly weird
%12--atnention very nice
% 13 -- early attn effect
% 14 -- attn, interesting interaction with saccade, few trials?
% 15 -attn strong, look closer at saccade
%16 --attn--weird interaction
%17-- good attn and saccade tuning
%18-- huge attn, huge saccade (example) %fix cluster issue, should be 1
%19--huge attn, huge saccade (example) 135, map mediocre
%20--ok, not clear, early attn
%21--ok, not clear
%22-- attn, not clear saccade


cell_file_dir='C:\Users\Sara\Desktop\DataAnalysis\cell_file_manager';
load(fullfile(cell_file_dir,'LIP_data.mat'));
{LIP_cell.dir};

cell_no=18; %18 vg, 12 ok.16 g
base_dir='C:\Users\Sara\Desktop\DataAnalysis\Michel';
load(fullfile(base_dir,LIP_cell(cell_no).dir,'proc',LIP_cell(cell_no).MGS_file.mat{1}));
ifname=fullfile(base_dir,LIP_cell(cell_no).dir,'proc',LIP_cell(cell_no).MGS_file.hdf{1});



bid = 1; % Channel number (1 for Blackrock recordings in Freiwald lab)
tres = 1e9/30000; % Specified in nanoseconds, reciprocal of sample rate
clussel =LIP_cell(cell_no).unit; % Cluster number(s)
make_grand_psth=1;
overwrite_grand_psth=1;

psth_mode='use_single_trials';%DO NOT CHANGE %'borisov_mean';% 'borisov_sem';


psthfilter.type = 'Gauss';
psthfilter.fs        = round(1e9/tres); %Sampling frequency (Hz)
psthfilter.sigma_ms  = 25; %Specified in milliseconds
psthfilter.length_ms = 100;   %Specified in milliseconds, truncation length, usually 4x more than sigma
psthdec = 100;


start_trig(1).type ='Marker';%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig(1).name ='Stimulus_onset';%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig(1).offset = -0.5; %-.05%1.0%-.25%-.05; %usually a negative value
end_trig(1).type ='Marker';%'TrialEnd'; %'Marker'%%
end_trig(1).name ='Fixation_out';%; %'Stimulus_onset';%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig(1).offset =-.5;%1.8;% 4.0%.05;%3.0


start_trig(2).type ='Marker';%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig(2).name ='Stimulus_onset';%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig(2).offset = -1.5; %-.05%1.0%-.25%-.05; %usually a negative value
end_trig(2).type ='Marker';%'TrialEnd'; %'Marker'%%
end_trig(2).name ='Fixation_out';%; %'Stimulus_onset';%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig(2).offset =.2;%1.8;% 4.0%.05;%3.0

% ------------ END OF USER-EDITABLE BLOCK -----------------
clear dhfun
% Design the filter
psthfilter = design_psth_filter(psthfilter);

fid = dhfun('open',ifname);
tmap = dh_get_trialmap_struct(fid);
%% Makes all relevant data structures
trls = [tds{1}.trials]; %holds all the basic trial parameters
%this gets you the name, number, rho and phi for the target surface of
%every trial
%[surf_str] = surface_info_withcatch(tmap,tds);
%[surf_str] = trial_info(tmap,tds);
%surface_numbers = [surf_str.numb];


%getting important info for MGS trials

targ_names = [tds{1}.doc_data(14).OBJECTPARAMS_SACMEM.objectName];

for j= 1:length(trls);
   mgs_str(j).tname =  find(strcmp(trls(j).object.name,targ_names)==1);
   mgs_str(j).t_idx =  find(strcmp(trls(j).object.name,targ_names)==1);
   mgs_str(j).tphi =  str2double(tds{1}.doc_data(j).OBJECTPARAMS_SACMEM(mgs_str(j).t_idx).posPhi);
   mgs_str(j).outcome = tmap.oc(j); 
end
%This is Michael Code (partially redundant):----

%td = tds{1};
%nocue = false;
%surf_str_extra = heiko_log_surf_params(td,nocue); % used to be called p in michael code
%---
%Sanity checks:
% if any(strcmp({surf_str.out},'Success')'~=(tmap.oc==7))
%     error('Tmap and surface info are not consistent')
% end
% test1=surf_str_extra.cuedsurfnames;
% test2=[surf_str.name];
% if ~all(strcmp(test1,test2))
%     error('The surf_str surface info structures are not consistent')
% end
%-----

%Load the spikes:

nspk = double(dhfun('GETNUMBERSPIKES',fid,bid)); %the numbers of spikes (total)
tstamps = dhfun('READSPIKEINDEX',fid,bid,1,nspk); %times that every spike occured

[spike_clus] = dhfun('READSPIKECLUSTER',fid,bid,1,nspk); %tells you cluster information
spks = tstamps(spike_clus == clussel);


%% Making Single-trial PSTHs for finer analysis.
if make_grand_psth
    if ~exist(fullfile(base_dir,LIP_cell(cell_no).dir,'proc',['cell_' sprintf('%03.g',cell_no) '_single_trial_mgs_PSTH.mat']),'file') | overwrite_grand_psth
        grand_psth_ingredients.fid=fid;
        grand_psth_ingredients.tmap=tmap;
        grand_psth_ingredients.str=mgs_str;
        grand_psth_ingredients.spks=spks;
        grand_psth_ingredients.tres=tres;
        grand_psth_ingredients.psthfilter=psthfilter;
        grand_psth_ingredients.psthdec=psthdec;
                
        for align=1:2 % Adapt the following code!!!!
            grand_psth_ingredients.align=align;
            grand_psth_ingredients.start_trig=start_trig(align);
            grand_psth_ingredients.end_trig=end_trig(align);
            
            [output_matrix, output_trials_used]= make_grand_psth_matrix_mgs(grand_psth_ingredients);
            grand_psth.matrix{align}=output_matrix;
            grand_psth.trials_used{align}=output_trials_used;
        end
        
        grand_psth.start_trig=start_trig;
        grand_psth.end_trig=end_trig;
        grand_psth.psthfilter=psthfilter;
        grand_psth.psthdec=psthdec;
        grand_psth.cell_no=cell_no;
        save(fullfile(base_dir,LIP_cell(cell_no).dir,'proc',['cell_' sprintf('%03.g',cell_no) '_single_trial_attn_PSTH.mat']),'grand_psth')
    else
        load(fullfile(base_dir,LIP_cell(cell_no).dir,'proc',['cell_' sprintf('%03.g',cell_no) '_single_trial_attn_PSTH.mat']),'grand_psth')
    end
end


trials_used=grand_psth.trials_used{1}& grand_psth.trials_used{2};
%% Make PSTHs of all interaction combinations
%Find all combs of phi positions and brt:
conditions_plotted=allcomb(unique([mgs_str.tphi]));
%Repeat the output with columns interchanged for plotting w/fixed brt
%instead of fixes phi:
conditions_plotted=[conditions_plotted; conditions_plotted(:,end:-1:1)];
%Number of unique phi (or brt) values:
cycle_length=length(unique([mgs_str.tphi]));


number_of_figures=2; %One for fixed attention, one for fixed saccade target
f_handles=zeros(number_of_figures,1);
linehandles=cell(cycle_length,2,2);%subplot_row,suplot_column,figure_number
ylims=zeros(size(f_handles));
legend_values=cellfun(@num2str,num2cell((unique([mgs_str.tphi])')),'UniformOutput',false);
for a=1:size(f_handles,1)
    for b=1:size(f_handles,2)
        f_handles(a,b)=figure;
        set(f_handles(a,b),'Position',get(0,'ScreenSize'));
    end
end
colorlist1=get(gca,'ColorOrder');
colorlist = [colorlist1; 0.75 0 0.75];

for j=1:length(conditions_plotted)
this_plot_cond=nan(length(trials_used),size(conditions_plotted,2));
this_plot_cond(:,1)= [mgs_str.tphi]==conditions_plotted(j,1);
%this_plot_cond(:,2)= [surf_str.brt]==conditions_plotted(j,2);

this_plot_trials=logical(trials_used.*prod(this_plot_cond,2));


for align=1:2
%     subplot_row=mod(ceil(j/(cycle_length))-1,cycle_length)+1;
     fig_in_use=gcf; %ceil(j/(length(conditions_plotted)/2));
%     sp_position=mod(ceil(j/(cycle_length))-1,cycle_length)*2+align;
     color_plotted=mod(j-1,cycle_length)+1;
%     set(0,'currentFigure',f_handles(fig_in_use));
%     %set(f_handles(ceil(j/(length(conditions_plotted)/2))),'Position',get(0,'ScreenSize'));
%     
%     
%     
%     subplot(cycle_length,2,sp_position)
%     
    switch psth_mode
        case 'borisov_sem' %obsolete
            if align==1
                ta_vect=ts_vect;
            elseif align==1
                ta_vect=te_vect;
            end
            warning('This psth_mode uses dh_calc_psth_stderr, which is slow and has yet to handle variable alignment points. PP 5/23/13')
            psth=calc_psth(spks,ts_vect(this_plot_trials),te_vect(this_plot_trials),ta_vect(this_plot_trials),tres,psthfilter.Num,1,psthfilter.delay,psthdec);
            t = 0:tres*psthdec/1e9:(length(psth.data)-1)*tres*psthdec/1e9; %time in seconds
            psth_sem=dh_calc_psth_stderr(fid,bid,ts_vect(this_plot_trials),te_vect(this_plot_trials),clussel,tres,psthfilter.Num,psthfilter.Den,psthfilter.delay,psthdec);
            %psth_se=dh_calc_psth_stderr(fid,bid,ts_vect(this_plot_trials),te_vect(this_plot_trials),ta_vect(this_plot_trials),clussel,tres,psthfilter.Num,psthfilter.Den,psthfilter.delay,psthdec);
            legendlist = ['test'];
            color2 = [1 0 0];
            [maxvalue minvalue]= shadowcaster_ver3(t', psth.data', 2*psth_sem.stderr', legendlist,color2);
            
        case 'borisov_mean' %obsolete
            if align==1
                ta_vect=ts_vect;
            elseif align==1
                ta_vect=te_vect;
            end
            psth=calc_psth(spks,ts_vect(this_plot_trials),te_vect(this_plot_trials),ta_vect(this_plot_trials),tres,psthfilter.Num,1,psthfilter.delay,psthdec);
            t = 0:tres*psthdec/1e9:(length(psth.data)-1)*tres*psthdec/1e9; %time in seconds
            psth_plot(psth,tres*psthdec,[],color2,[1 0 0],[],bid,clussel,'Spikes per second');
            %psth_plot(psth,tres,dopts,gcolor,contrib_color,ax,elec,clust,ylab)
            
            %Line indicating surface onset (stmulus onset+pause_duration):
            line([-start_trig.offset+unique(surf_str_extra.pausedur(trials_used)) -start_trig.offset+unique(surf_str_extra.pausedur(trials_used))], ylim, 'Color', 'k', 'LineStyle', '--')
            %Line indicating cue onset (the actual stimulus onset):
            line([-start_trig.offset -start_trig.offset], ylim, 'Color', [.3 .3 .3], 'LineStyle', '--')
            %Make xlim a reasonable range: (0.5 sec before cue onset until 2sec after surface onset)
            xlim_used=[-start_trig.offset-0.5 -start_trig.offset+unique(surf_str_extra.pausedur(trials_used))+2 ];
            xlim(xlim_used)
            
        case 'use_single_trials' %current approach
            plot_params.current_yrange=ylims(fig_in_use);
            plot_params.align=align;
            plot_params.tres=tres;
            plot_params.legendlist=['test'];           
            plot_params.color=colorlist(color_plotted,:);
            [linehandle,current_ylim]=make_grand_psth_plot_mgs(grand_psth,this_plot_trials,plot_params);
            %linehandles{subplot_row,align,fig_in_use}=[linehandles{subplot_row,align,fig_in_use}; linehandle];
            ylims(fig_in_use)=current_ylim;
        otherwise
            error('invalid psth mode')
    end
    if fig_in_use==1
      %  title([' Attention at phi=' num2str(conditions_plotted(j,1))]);
    else
%        title(['Saccade at phi=' num2str(conditions_plotted(j,2))]);

    end
   % lg=legend(linehandles{subplot_row,align,fig_in_use},legend_values(1:color_plotted),'Location','NorthWest');
%     set(lg,'FontSize',5.5,'Box','off')
%     set(lg,'units','pixels'); 
%     lp=get(lg,'outerposition');
%     new_width=20; new_height=35;
%     set(lg,'outerposition',[lp(1),lp(2)+0.9*(lp(4)-new_height),new_width,new_height]);



    
end
end

titles={'Saccade Effect: One attention condition per row','Attention Effect: One saccade target per row'}
for i=1:length(f_handles)
set(0,'currentFigure',f_handles(i));
[axh,labelh]=suplabel(titles{i},'t');
set(labelh,'FontSize',15)
end

%% Makes plot for main effects:

%Main effect of attention:
Main_effects(1).name='Target';
Main_effects(1).values=unique([mgs_str.tphi]);
Main_effects(1).variable_used='[mgs_str.tphi]';

Main_effects(2).name='Saccade';
Main_effects(2).values=unique([mgs_str.tphi]);
Main_effects(2).variable_used='[mgs_str.tphi]';
titles={'Aligned on Stimulus Onset','Aligned on Saccade Onset'}
for i=1:length(Main_effects)
   f=figure;
   set(f,'Position',get(0,'ScreenSize'));
   eval(['var_used=' Main_effects(i).variable_used ''';']);
   linehandles=cell(1,2);
   plot_params.current_yrange=0;
   for j=1:length(Main_effects(i).values)       
       for align=1:2
           subplot(1,2,align);
           title(titles{align})
           hold on
           this_plot_trials= trials_used & (var_used==Main_effects(i).values(j));
           plot_params.align=align;
           plot_params.tres=tres;
           plot_params.legendlist=[];
           plot_params.color=colorlist(j,:);
           [linehandle,curent_yarange]=make_grand_psth_plot_mgs(grand_psth,this_plot_trials,plot_params);
           plot_params.current_yrange=curent_yarange;
           linehandles{align}=[linehandles{align};linehandle];
           if j==length(Main_effects(i).values)
               lg=legend(linehandles{align},legend_values,'Location','SouthWest');
               set(lg,'Box','off')
                   set(lg,'units','pixels'); 
    lp=get(lg,'outerposition');
    set(lg,'outerposition',[lp(1:2),30,lp(4)]);
           end
           
       end
   end
   [axh,labelh]=suplabel(['Main effect of ' Main_effects(i).name],'t');
   set(labelh,'FontSize',15)
end
% 
% %% Make plot by physcial conditions
% locations=unique([surf_str.phi]);
% locations=reshape(locations,2,2); %rows indicate conditions
% 
% f=figure;
% set(f,'Position',get(0,'ScreenSize'));
% 
% colorlist=get(gca,'ColorOrder');
% colorlist=colorlist(1:3,:);
% colorlist(end+1,:)=0.25*[1 1 1];
% ylims=0;
% for surf_cond=1:2
%     surf_pos=locations(surf_cond,:)
%     for sacc_cond=1:2
%         sacc_pos=locations(sacc_cond,:)
%         for align=1:2
%             subplot_row_used= 2*(surf_cond-1) +sacc_cond;
%             set(0,'currentFigure',f);
%             subplot_used=(subplot_row_used-1)*2+align;
%             subplot(4,2,subplot_used)
%             linehandles=[];
%             legendtexts={};
%             for surf_p_used=1:length(surf_pos)
%                 for sacc_p_used=1:length(sacc_pos)
%                     this_plot_trials=trials_used & ([surf_str.phi]'==surf_pos(surf_p_used)) ...
%                         & ([surf_str.brt]'==sacc_pos(sacc_p_used));
%                     color_plotted=2*(surf_p_used-1)+sacc_p_used;
%                     
%                     plot_params.current_yrange=ylims;               
%                     plot_params.align=align;
%                     plot_params.tres=tres;
%                     plot_params.legendlist=[];
%                     plot_params.color=colorlist(color_plotted,:);
%                     [linehandle,current_ylim]=make_grand_psth_plot(grand_psth,this_plot_trials,surf_str_extra,plot_params);
%                     linehandles(end+1)=linehandle;
%                     legendtexts{end+1}=['s=' num2str(surf_pos(surf_p_used)) '; t=' num2str(sacc_pos(sacc_p_used))];
%                     ylims=current_ylim;
%                     
%                 end
%             end
%             lg=legend(linehandles,legendtexts,'Location','NorthWest');
%             set(lg,'FontSize',7,'Box','off')
%             set(lg,'units','pixels');
%             lp=get(lg,'outerposition');
%             new_width=30; new_height=35;
%             set(lg,'outerposition',[lp(1),lp(2)+0.9*(lp(4)-new_height),new_width,new_height]);
%             
%             
%         end
%     end
% end
%% Make tuning curves



%% Close HDF file  
dhfun('close',fid);


