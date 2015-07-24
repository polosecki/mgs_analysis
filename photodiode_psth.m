function [grand_psth_matrix, trials_used, time_axis] = make_grand_psth_matrix_mgs(grand_psth_ingredients)
%Makes a matrix with individual trial PSTHs:
% Required fields are:
% fid
% tmap
% start_trig
% end_trig
% str %structure with trial info
%align %=1 for stim onset or 2 for saccade onset
%spks
%tres
%psthfilter

fid=grand_psth_ingredients.fid;
tmap=grand_psth_ingredients.tmap;
start_trig=grand_psth_ingredients.start_trig;
end_trig=grand_psth_ingredients.end_trig;
outcomes = [grand_psth_ingredients.str.outcome]';
spks=grand_psth_ingredients.spks;
tres=grand_psth_ingredients.tres;
psthfilter=grand_psth_ingredients.psthfilter;
psthdec=grand_psth_ingredients.psthdec;
align=grand_psth_ingredients.align;

%Get time start and end (in nanoseconds) for each trial:
[time_start,time_end,inclidx] = find_ref_points(fid,tmap,start_trig,end_trig);

%Prevent future errors!
if find(time_start<0)
    warning('Hey, your time_start is less than zero!')
end

trials_with_good_markers=zeros(length(outcomes),1);
trials_with_good_markers(inclidx)=1;
ts_vect=nan(length(trials_with_good_markers),1);
te_vect=nan(length(trials_with_good_markers),1);
ts_vect(logical(trials_with_good_markers))=time_start;
te_vect(logical(trials_with_good_markers))=time_end;
durations=te_vect-ts_vect;
trials_used = zeros(length(outcomes),1);
trials_used = (trials_with_good_markers==1);

%Use max(durations) to predefine the size of the nan matrix where psth will
%occur.
[~,trialID]=max(durations.*trials_used); %the max durations of the trials used

longest_psth = calc_psth(spks,ts_vect(trialID),te_vect(trialID),ts_vect(trialID),tres,psthfilter.Num,1,psthfilter.delay,psthdec);

grand_psth_matrix=nan(length(trials_used),length(longest_psth.data));
initial_size=size(grand_psth_matrix,2);
for trialID=1:length(trials_used)
    if trials_used(trialID)
        psth=calc_psth(spks,ts_vect(trialID),te_vect(trialID),ts_vect(trialID),tres,psthfilter.Num,1,psthfilter.delay,psthdec);
        if align==1
            grand_psth_matrix(trialID,1:length(psth.data))=psth.data;
        elseif align==2
            grand_psth_matrix(trialID,end-(length(psth.data)-1):end)=psth.data;
        else
            error('Alignment mode not supported')
        end
    end
end
if size(grand_psth_matrix,2)~=initial_size
    error('Possible trial missalignment');
end

time_axis = 0:tres*psthdec/1e9:(initial_size-1)*tres*psthdec/1e9; %time in seconds
if  align==1
    %Center on surface onset(stmulus onset+pause_duration):
    %t_zero=-grand_psth.start_trig(plot_params.align).offset+unique(surf_str_extra.pausedur(this_plot_trials));
    t_zero=-start_trig.offset; %for mgs, no pause duration
elseif align==2
    %Center on saccade onset:
    t_zero=time_axis(end)-end_trig.offset;
end
    time_axis=time_axis-t_zero;
