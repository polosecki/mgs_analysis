function [linehandles]= shadowcaster_ver3(x, meanvalue, sem, legendlist,colorlist)
%plot mean values along provided x label, with error shades
%no additional function required
%ver3: Akinori F. Ebihara, June 29, 2011.
%ver3PP: Modified by PP (05/24/2013) to provide a handle to the plotted lines as output .

edge = 'none'; %edge color
add=0; %add=1 is to add to current plot. otherwise set 0
transparency=0.15;

minvalue = Inf; maxvalue = -Inf;
linehandles = [];
for i = 1:size(meanvalue,2)
    UpperLimit = meanvalue(:,i)+sem(:,i);
    LowerLimit = meanvalue(:,i)-sem(:,i);
    
    if minvalue > min(min(LowerLimit))
        minvalue = min(min(LowerLimit));
    end
    if maxvalue < max(max(UpperLimit))
        maxvalue = max(max(UpperLimit));
    end
    
    color = colorlist(i,:);
    
if length(LowerLimit)==length(x)
    ypoints=[UpperLimit;flipud(LowerLimit)]; %orignal from ae
    xpoints = [x;flipud(x)];
    if add
        hold on
    end
    fillhandle=fill(xpoints,ypoints,color);
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);
    if add
        hold off
    end
else
    sprintf('Error: inconsistent argument length')
end
    hold on
    hand = plot(x,meanvalue(:,i),'Color',color);
    linehandles = [linehandles, hand];
end
if nargin==4
    legend(linehandles,legendlist,'Location','NE');
end