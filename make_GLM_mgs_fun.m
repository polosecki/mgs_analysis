function [results]=make_GLM_mgs_fun(cell_no,monkey,area,operation_mode,z_score_data,cross_time_comparisons)

%Inputs: cell_no: cell number in data base
%         monkey: monkey name
%         area: brain area
%         operation_mode: 'time_course' for full time analysis with plot
%                         'fixed_points' for just fixed set of time points (for histogram of population)
%         z_score_data: 

%cell_no=28;%PITd; %28;% LIP;


base_dir='/Freiwald/ppolosecki/lspace';
%if strcmp(monkey,'Michel') & strcmp(area,'LIP')
%    base_dir='/Freiwald/ppolosecki/lspace/sara';
%end
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
%monkey='Quincy';
%area='LIP';%'IP';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

fixed_time_bins={[0.15]; %stim-onset;
                 [-0.5 -.1]; %one seclong before saccade onset
                 }; % 

%% Load basic files

load(cell_file);
load(results_file)
load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).MGS_file.mat{end}));
ifname=fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).MGS_file.hdf{end});

fid = dhfun('open',ifname);
trls = [tds{1}.trials];
tmap = dh_get_trialmap_struct(fid);

targ_names = [tds{1}.doc_data(1).OBJECTPARAMS_SACMEM.objectName];

for j= 1:length(trls);
   mgs_str(j).tname =  find(strcmp(trls(j).object.name,targ_names)==1);
   mgs_str(j).t_idx =  find(strcmp(trls(j).object.name,targ_names)==1);%what is the diffrence with above line??
   mgs_str(j).tphi =  str2double(tds{1}.doc_data(j).OBJECTPARAMS_SACMEM(mgs_str(j).t_idx).posPhi);
   mgs_str(j).outcome = tmap.oc(j); 
end

%---

%-----

if ~exist(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_mgs_PSTH.mat']),'file');
    error('Run the first atention analysis code to create the PSTHs')
else
    load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_mgs_PSTH.mat']));
end


%% Run GLM and Tests:
%Commands of possible relevance:
%
%mvregress: multivariate general linear model
%manova1:multavariate analysis of variance
%fast_glmfit: fit glm model (FSFAST)
%fast_fratio: evaluate contrasts
%anovan: N-way anova


num_matrices=2;


if z_score_data
    good_times={[-0.5 0.5],[-.8 0.2]};
    raw_grand=[];
    for i=1:num_matrices
        idx=(grand_psth.time_axis{i}>good_times{i}(1) & grand_psth.time_axis{i}<good_times{i}(2));
        raw_grand=[raw_grand grand_psth.matrix{i}(:,idx)];
    end
    conditions=unique([mgs_str.tphi]');
    mean_mat=nan(size(conditions,1),size(raw_grand,2));
    zz=zeros(length(grand_psth.trials_used{1}),num_matrices);
    for col=1:num_matrices
        zz(:,col)=grand_psth.trials_used{col};
    end
    trials_used=prod(zz,2);
    for j=1:length(conditions)
        these_trials=trials_used & [mgs_str.tphi]'==conditions(j);
        mean_mat(j,:)=nanmean(raw_grand(these_trials,:));
    end
    
    mean_center=nanmean(mean_mat(:));
    std_scale=nanstd(mean_mat(:));
else
    mean_center=0;
    std_scale=1;
end
 
for mat_used=1:num_matrices
    

    t=grand_psth.time_axis{mat_used};
    %trial_count=sum(~isnan(grand_psth.matrix{mat_used}),1);
    %times_used=trial_count>15;
    RF_surf=nan(length(results_data),1);
    RF_surf(~cellfun(@isempty,{results_data.closest_surf_phi}))=[results_data.closest_surf_phi]';
    
    switch operation_mode
        case 'time_course'
            tbins_center=t(3:3:end); %50ms bins
            tbins_semi_width=0.075; % i.e., width=2*semi_width
        case 'fixed_points'
            tbins_center=fixed_time_bins{mat_used};
            tbins_semi_width=0.1; %200 ms width
        case 'betas_for_pca'
            tbins_center=t;
            tbins_semi_width=0;
        otherwise
            error('Set a valid Operation Mode')
            
    end
    y=nan(size(grand_psth.matrix{mat_used},1),length(tbins_center));
    for i=1:length(tbins_center)
        y(:,i)=nanmean(grand_psth.matrix{mat_used}(:,abs(t-tbins_center(i))<=tbins_semi_width),2);
    end

     y=(y-mean_center)/std_scale;

    [temp]= make_GLM_and_contrasts_from_inst_firing_mgs(y,RF_surf(cell_no),[mgs_str.tphi]');
    switch operation_mode
        case 'fixed_points'
            if cross_time_comparisons & mat_used==2
                cross_time_results=make_cross_time_GLM_mgs(y,RF_surf(cell_no),[mgs_str.tphi]');
                temp=[temp cross_time_results];
            end
    end
    
    results{mat_used}=temp; clear temp;
    results{mat_used}(1).time=tbins_center;
    results{mat_used}(1).y=y;
end
%% Make plots
switch operation_mode
    case 'time_course'
        contrasts_plotted={logical([1])};
        plot_GLM_contrasts_mgs(results,contrasts_plotted,cell_str,cell_no);
    case 'betas_for_pca'
        contrasts_plotted={logical([0 0 0 0 0 0 0 0 0 0]);
                           logical([0 0 0 0 0 0 0 0 0]);
                           logical([1 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0])};
                       %plot_GLM_contrasts(results,contrasts_plotted);
end


