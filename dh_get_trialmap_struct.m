function [tmap] = dh_get_trialmap_struct(fid)

[tn,sn,oc,ts,te] = dhfun('GETTRIALMAP',fid);

tn = double(tn);
sn = double(sn);
oc = double(oc);

tmap.tn = tn;
tmap.sn = sn;
tmap.oc = oc;
tmap.ts = ts;
tmap.te = te;
% Index of trial numbers - will be updated during trial selection.
% Index is used for modifications of the trialmap - to restore
% the original position of a certain trial in the whole trialmap.
tmap.idx = 1:length(tn);
