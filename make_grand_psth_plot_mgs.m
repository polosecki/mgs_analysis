function [linehandle,curent_yarange]=make_grand_psth_plot_mgs(grand_psth,this_plot_trials,plot_params)

grand_matrix=grand_psth.matrix{plot_params.align};


psth=nanmean(grand_matrix(this_plot_trials,:));
psth_std=nanstd(grand_matrix(this_plot_trials,:));
contrib=sum(~isnan(grand_matrix(this_plot_trials,:)));
use_idx=contrib>1;
psth_sem=psth_std;
psth_sem(use_idx)=psth_std(use_idx)./sqrt(contrib(use_idx));

tres=plot_params.tres;
psthdec=grand_psth.psthdec;
t = 0:tres*psthdec/1e9:(length(psth)-1)*tres*psthdec/1e9; %time in seconds
if plot_params.align==1
    %Center on surface onset(stmulus onset+pause_duration):
    %t_zero=-grand_psth.start_trig(plot_params.align).offset+unique(surf_str_extra.pausedur(this_plot_trials));
    t_zero=-grand_psth.start_trig(plot_params.align).offset; %for mgs, no pause duration
elseif plot_params.align==2
    %Center on saccade onset:
    t_zero=t(end)-grand_psth.end_trig(plot_params.align).offset;
end
    t=t-t_zero;
use_idx=contrib>0;
[linehandle]= shadowcaster_ver3PP(t(use_idx)', psth(use_idx)', 2*psth_sem(use_idx)', plot_params.legendlist,plot_params.color);

if plot_params.align==1
    %Line indicating surface onset (stmulus onset+pause_duration):
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--')
    %Line indicating cue onset (the actual stimulus onset):
    %line([-unique(surf_str_extra.pausedur(this_plot_trials)) -unique(surf_str_extra.pausedur(this_plot_trials))], ylim, 'Color', [.3 .3 .3], 'LineStyle', '--')
    %Make xlim a reasonable range: (0.5 sec before cue onset until 2sec after surface onset)
    xlim_used=[-0.5 0.8];
elseif plot_params.align==2
    %Line indicating saccade onset:
%     line([t(end)-grand_psth.end_trig.offset t(end)-grand_psth.end_trig.offset], ylim, 'Color', 'k', 'LineStyle', '--')
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--')
    %Make xlim a reasonable range: (2 sec before saccade until 0.5 sec after saccade)
    xlim_used=[-0.8 0.5];
end
xlim(xlim_used)
range_y(1)=0;
[C,bin_start]=min(abs(t-xlim_used(1)));
[C,bin_end]=min(abs(t-xlim_used(2)));
range_y(2)=max(max(psth(bin_start:bin_end)+2*psth_sem(bin_start:bin_end)),plot_params.current_yrange);
curent_yarange=range_y(2);
ylim(range_y)
box off