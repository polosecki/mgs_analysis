function [linehandles]=plot_GLM_contrasts(results,contrasts_plotted,cell_str,cell_no)
%INPUTS:
%   results: results structure as provided by make_GLM_and_contrasts_from_inst_firing
%   contrasts_plotted: 2-by-1 cell with logical vectors indicating the contrasts to be plotted from


    num_matrices=2;
    subplot_column_used=[1 2];


show_legends=1;
f=figure;
set(f,'Position',get(0,'ScreenSize'));
lims_used=zeros(num_matrices,2);
%colors_used=distinguishable_colors(sum(contrasts_plotted{1}+sum(contrasts_plotted{2})));
colors_used=distinguishable_colors(sum(contrasts_plotted{1}));

for mat_used=1:num_matrices
    subplot(1,num_matrices,subplot_column_used(mat_used)); hold on
    if mat_used==1
        xlim([-0.5 0.8])
    elseif mat_used==2
        xlim([-0.8 0.025]);%0.5])
    elseif mat_used==3
    else
        error('what mat??')
    end
    
    
    y_matrix=[results{mat_used}.GLM(1).ces(contrasts_plotted{1},:)];
    y__std_matrix=[results{mat_used}.GLM(1).ces_std(contrasts_plotted{1},:)];
    [temp]= shadowcaster_ver3PP(results{mat_used}.time',y_matrix', 2*y__std_matrix', [],colors_used);
    if mat_used==1
        linehandles=temp;
    end
    clear temp;
    lims_used(mat_used,:)=ylim;
end

for mat_used=1:num_matrices
    subplot(1,num_matrices,subplot_column_used(mat_used))
    ylim([min(lims_used(:,1)) max(lims_used(:,2))]);
    line(xlim,[0 0],'Color','black','LineStyle','-')
    line([0 0],ylim,'Color','black','LineStyle','--')
    xlabel ('Time from trigger event(s)')
    if mat_used==1
        ylabel('Effect Size (z-score units)')
        if show_legends
%             a={results{1}.GLM(1).contrast.name{contrasts_plotted{1}} results{1}.GLM(2).contrast.name{contrasts_plotted{2}}};
%             for hh=1:length(a)
%                 disp([num2str(hh) ' ' a{hh}]);
%             end
        lg=legend(linehandles,results{1}.GLM(1).contrast.name{contrasts_plotted{1}});
        set(lg,'Location','SouthWest');
        set(lg,'Box','off')
        %         set(lg,'units','pixels');
        %         lp=get(lg,'outerposition');
        %         set(lg,'outerposition',[lp(1:2),50,lp(4)]);
        end
    end
end

if nargin>3
[axh,labelh]=suplabel([cell_str(cell_no).monkey '__' cell_str(cell_no).area '__cell_no:__' sprintf('%02.0f',cell_no)],'t');
set(labelh,'Interpreter','none');
end

