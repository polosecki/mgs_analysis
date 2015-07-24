
% Usage:
%  [] = mgs_analysis_func_pp(cell_file,cell_no,basedir,overwrite_grand_psth,figure_dir)
% Inputs:
%  cell_file: the file with the strcuture (with its path)
%  cell_no: the cell number
%  basedir: the directory where actuall datafiles are stored
%  overwrite_grand_psth: 1 if you want to overwrite rxisting grand_psth file
%  figure_dir: directory to save figures (leave empty if dont want to save)
% Output:
%  none
close all; clear all;

figure_dir=[];
closest_surf_phi=0;
%load(fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).MGS_file.mat{end}));
%ifname=fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).MGS_file.hdf{end});
load('/Freiwald/ppolosecki/lspace/quincy/pablo_mgs_length/proc/pablo_measurements_0000.mat');
ifname='/Freiwald/ppolosecki/lspace/quincy/pablo_mgs_length/proc/polosecki_mgs_duration_001_r.dh5'


bid = 137; % Channel number (1 for Blackrock recordings in Freiwald lab)
tres = 1e9/30000; % Specified in nanoseconds, reciprocal of sample rate
clussel =0; % Cluster number(s)
make_grand_psth=1;
overwrite_grand_psth=1;
save_figures=~isempty(figure_dir);


psthfilter.type = 'Gauss';
psthfilter.fs        = round(1e9/tres); %Sampling frequency (Hz)
psthfilter.sigma_ms  = 1; %Specified in milliseconds
psthfilter.length_ms = 4;   %Specified in milliseconds, truncation length, usually 4x more than sigma
psthdec = 100;


start_trig(1).type ='Marker';%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig(1).name ='Stimulus_onset';%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig(1).offset = -0.5; %-.05%1.0%-.25%-.05; %usually a negative value
end_trig(1).type ='Marker';%'TrialEnd'; %'Marker'%%
end_trig(1).name ='Stimulus_onset';%; %'Stimulus_onset';%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig(1).offset =.8;%1.8;% 4.0%.05;%3.0


start_trig(2).type ='Marker';%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig(2).name ='Stimulus_onset';%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig(2).offset = 0.4; %-.05%1.0%-.25%-.05; %usually a negative value
end_trig(2).type ='Marker';%'TrialEnd'; %'Marker'%%
end_trig(2).name ='Fixation_out';%; %'Stimulus_onset';%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig(2).offset =.2;%1.8;% 4.0%.05;%3.0

% ------------ END OF USER-EDITABLE BLOCK -----------------
addpath(fullfile('.','plot2svg'))
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

targ_names = [tds{1}.doc_data(1).OBJECTPARAMS_SACMEM.objectName];

for j= 1:length(trls);
   mgs_str(j).tname =  find(strcmp(trls(j).object.name,targ_names)==1);
   mgs_str(j).t_idx =  find(strcmp(trls(j).object.name,targ_names)==1);%what is the diffrence with above line??
   mgs_str(j).tphi =  str2double(tds{1}.doc_data(j).OBJECTPARAMS_SACMEM(mgs_str(j).t_idx).posPhi);
   mgs_str(j).outcome = tmap.oc(j); 
end


%Load the spikes:

nspk = double(dhfun('GETNUMBERSPIKES',fid,bid)); %the numbers of spikes (total)
tstamps = dhfun('READSPIKEINDEX',fid,bid,1,nspk); %times that every spike occured

[spike_clus] = dhfun('READSPIKECLUSTER',fid,bid,1,nspk); %tells you cluster information
spks = tstamps(spike_clus == clussel);


%% Making Single-trial PSTHs for finer analysis.
if make_grand_psth
    if  overwrite_grand_psth
        grand_psth_ingredients.fid=fid;
        grand_psth_ingredients.tmap=tmap;
        grand_psth_ingredients.str=mgs_str;
        grand_psth_ingredients.spks=spks;
        grand_psth_ingredients.tres=tres;
        grand_psth_ingredients.psthfilter=psthfilter;
        grand_psth_ingredients.psthdec=psthdec;
                
        for align=1:1
            grand_psth_ingredients.align=align;
            grand_psth_ingredients.start_trig=start_trig(align);
            grand_psth_ingredients.end_trig=end_trig(align);
            
            [output_matrix, output_trials_used, output_time_axis]= photodiode_psth(grand_psth_ingredients);
            grand_psth.matrix{align}=output_matrix;
            grand_psth.trials_used{align}=output_trials_used;
            grand_psth.time_axis{align}=output_time_axis;
        end
        
        grand_psth.start_trig=start_trig;
        grand_psth.end_trig=end_trig;
        grand_psth.psthfilter=psthfilter;
        grand_psth.psthdec=psthdec;
        save('/Freiwald/ppolosecki/lspace/quincy/pablo_mgs_length/proc/matrix.mat','grand_psth')
    else
        load('/Freiwald/ppolosecki/lspace/quincy/pablo_mgs_length/proc/matrix.mat')
    end
end


trials_used=grand_psth.trials_used{1};

%% Center target position on RF-surface-centric coordinates
if isempty(closest_surf_phi)
    closest_surf_phi=180; %assume coutralateral RF;
end
angles=unique([mgs_str.tphi]);
RF_angles_distance=abs(angles-closest_surf_phi);
RF_angles_distance(RF_angles_distance>180)=RF_angles_distance(RF_angles_distance>180)-180;
[closest_targ_phi,min_index]=min(RF_angles_distance); 



%% Makes plot for main effects:
%Main effect ofsaacade target:
Main_effects(1).name='Target';
Main_effects(1).values=circshift(unique([mgs_str.tphi])',-min_index+1);
Main_effects(1).variable_used='[mgs_str.tphi]';

colorlist=distinguishable_colors(length(unique([mgs_str.tphi])));
%colorlist1=colorlist(1:4,:);colorlist2=colorlist(5:end,:);
%colorlist = reshape([colorlist1(:) colorlist2(:)]',2*size(colorlist1,1), []);
legend_values=cellfun(@num2str,num2cell((unique([mgs_str.tphi])')),'UniformOutput',false);
 
titles={'Aligned on Stimulus Onset','Aligned on Saccade Onset'};
for i=1:length(Main_effects)
   f=figure;
   set(f,'Position',get(0,'ScreenSize'));
   eval(['var_used=' Main_effects(i).variable_used ''';']);
   linehandles=cell(1,2);
   plot_params.current_yrange=0;
   for j=1:length(Main_effects(i).values)       
       for align=1:1
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
               lg=legend(linehandles{align},legend_values,'Location','NorthWest');
               set(lg,'Box','off')
                   set(lg,'units','pixels'); 
    lp=get(lg,'outerposition');
    set(lg,'outerposition',[lp(1:2),30,lp(4)]);
           end
           
       end
   end
      for z=1:2
        subplot(1,2,z)
        ylim([0 curent_yarange])
      end
   [axh,labelh]=suplabel(['Effect of ' Main_effects(i).name],'t');
   set(labelh,'FontSize',15)
%      [~,temp,~]=fileparts(cell_str(cell_no).MGS_file.hdf{end});
%    [axh,labelh]=suplabel([temp '____' cell_str(cell_no).monkey '__' cell_str(cell_no).area '__cell_no:__' sprintf('%02.0f',cell_no)],'x'); 
%    [axh,labelh]=suplabel(temp,'x');
%    set(labelh,'Interpreter','none');
    if save_figures
        filename=['cell_' sprintf('%03.0f',cell_no) '_mgs'];
        %saveas(f,fullfile(figure_dir,[filename '.eps']));
        plot2svg(fullfile(figure_dir,[filename '.svg']),f);
        saveas(f,fullfile(figure_dir,[filename '.fig']));
    end

end



%% Plot continuous photdiode data

elect=137;
sample_start=1;
sample_end=10;
ch_start=1;ch_end=1;

%dhfun('READCONT',fid,elect,sample_start,sample_end,ch_start,ch_start);


[time, offset]=dhfun('READCONTINDEX',fid,elect,1,1)
[PERIOD]=dhfun('GETCONTSAMPLEPERIOD',fid,elect)

[time_start,time_end,inclidx] = find_ref_points(fid,tmap,start_trig(1),end_trig(1));
start_indx=int32(round((time_start-double(time))/double(PERIOD)+1));

duration=int32(round(unique(time_end-time_start)/double(PERIOD)));

lfp_matrix=nan(length(tmap.ts),duration+1);

for i=1:size(lfp_matrix,1)
    lfp_matrix(i,:)=dhfun('READCONT',fid,elect,start_indx(i),start_indx(i)+duration,ch_start,ch_start);
end

t=start_trig(1).offset+double(PERIOD)/1e9*(0:(size(lfp_matrix,2)-1));
figure
imagesc(lfp_matrix)
figure;
plot(t,lfp_matrix')
figure;
plot(t,mean(lfp_matrix))

%% Close HDF file
dhfun('close',fid);


