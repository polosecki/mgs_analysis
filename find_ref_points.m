function [ts,te,inclidx] = find_ref_points(fid,tmap,start_trig,end_trig)
% Find the intervals of interest based on trigger points and offsets

global DH

excl = zeros(length(tmap.ts),1);

if strcmp(start_trig.type, 'TrialStart')
    % Just use the trialstart time which is present for all trials
    ts = tmap.ts;

elseif strcmp(start_trig.type,'TrialEnd')
    % Just use the trialend time which is present for all trials
    ts = tmap.te;

elseif strcmp(start_trig.type,'IntervalBegin')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,start_trig.name);
    [ts,excl] = proc_interval(tmap,its,ite,its,excl);
elseif strcmp(start_trig.type, 'IntervalEnd')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,start_trig.name);
    [ts,excl] = proc_interval(tmap,its,ite,ite,excl);
elseif strcmp(start_trig.type,'Marker')
    % Find markers which are contained by trials
    mt = dhfun('GETMARKER',fid,start_trig.name);
    [ts,excl] = proc_marker(tmap,mt,excl);
else
    error(['Unknown start_trig.type specification' start_trig.type]);
end

if strcmp(end_trig.type, 'TrialStart')
    te = tmap.ts;
    
elseif strcmp(end_trig.type, 'TrialEnd')
    te = tmap.te;
    
elseif strcmp(end_trig.type,'IntervalBegin')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,end_trig.name);
    [te,excl] = proc_interval(tmap,its,ite,its,excl);
elseif strcmp(end_trig.type,'IntervalEnd')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,end_trig.name);
    [te,excl] = proc_interval(tmap,its,ite,ite,excl);
elseif strcmp(end_trig.type,'Marker')
    % Find markers which are contained by trials
    mt = dhfun('GETMARKER',fid,end_trig.name);
    [te,excl] = proc_marker(tmap,mt,excl);
else
    error(['Unknown end_trig.type specification: ' end_trig.type]);
end

% Apply the offsets
ts = ts + round(1e9*start_trig.offset);
te = te + round(1e9*end_trig.offset);

% Exclude also the trials where te < ts
excl(find(~excl & (te < ts))) = 1;

% Now, exclude all trial intervals which were marked for exclusion
if any(excl)
    inclidx = find(~excl);
    ts = ts(inclidx);
    te = te(inclidx);
else
    inclidx = 1:length(ts);
end
% End of main function
% --------------------

function [tp,excl] = proc_interval(tmap,its,ite,itp,excl)
tp = zeros(length(tmap.ts),1);
for i=1:length(tmap.ts)
    idx = find(ite >= tmap.ts(i) & its <= tmap.te(i));
    % If there are no intervals found within this trial, or
    % more than one interval, exclude this trial
    if length(idx) ~= 1
        excl(i) = 1;
    else
        tp(i) = itp(idx);
    end
end

function [tp,excl] = proc_marker(tmap,mt,excl)
tp = zeros(length(tmap.ts),1);
for i=1:length(tmap.ts)
    idx = find(tmap.ts(i) <= mt & tmap.te(i) >= mt);
    % If there are no markers found within the trial, or
    % more than one marker, exclude this trial
    if length(idx) ~= 1
        excl(i) = 1;
    else
        tp(i) = mt(idx);
    end
end
