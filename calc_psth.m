function [psth] = calc_psth(tstamps,ts,te,ta,tres,Num,Den,fdelay,dec)
% Calculate the PSTH without loading data. Arguments are
% tstamps - spike time stamps
% ts, te - time intervals
% ta - alignment point times for each interval
% clussel - include only spikes with specified cluster number(s)
%           -1 means no cluster selection (all spikes)
% tres - time resolution (in nanoseconds)
% Num,Den - coefficients of FIR or IIR filter used
%           for histogram-like filtering
% fdelay - FIR filter delay to compensate for
%          (use zero for non linear-phase filters)
% dec - downsample the filtered and non-filtered result
%       by this factor

% Find the maximum time between trial start and alignment point
maxals = max(ta-ts);
% Find the offset of the alignment bin
align_offset = round(maxals/tres)+1;

% Find the maximum time between alignment point and trial end
maxale = max(te-ta);

% Find the number of bins in this trial
maxbw = align_offset + round(maxale/tres);

% Prepare the resulting psth array and supplemental data
contrib = zeros(1,maxbw);
psthdata = zeros(length(ts),maxbw);
bins = 1:maxbw;
%fprintf('maxbw = %d\n',maxbw);

% Process, trial by trial
for tk = 1:length(ts)
    odata = zeros(1,maxbw);
    otstamps = tstamps(tstamps>=ts(tk) & tstamps <= te(tk));
    % Refer all spikes to the alignment point (timestamps may get negative)
    oidx = round((otstamps-ta(tk))/tres)+align_offset;
    odata = hist(oidx,bins);

    bws = round((ts(tk)-ta(tk))/tres) + align_offset;
    bwe = round((te(tk)-ta(tk))/tres) + align_offset;
    contrib(bws:bwe) = contrib(bws:bwe) + 1;
    %fprintf('tk= %d tlenbw = %d lodata=%d\n',tk,tlenbw,length(odata));
    psthdata(tk,:) = odata;
end
% Apply the contribution info (where it is nonzero)
psthdata_nco = sum(psthdata,1);
idxcontr = find(contrib~=0);
stderrs = zeros(1,maxbw);
stdevs = zeros(1,maxbw);
stdevs(idxcontr) = sqrt(var(psthdata(:,idxcontr)));
stderrs(idxcontr) = sqrt(var(psthdata(:,idxcontr))./contrib(idxcontr));
psthdata = psthdata_nco(idxcontr) ./ contrib(idxcontr);
psthdata(contrib==0) = 0;

% Divide the values by time resolution, thereby
% obtaining the actual density of spikes per second
% (Firing rate)
psthdata = psthdata ./ tres .* 1e9; % tres is specified in nanoseconds
stderrs = stderrs ./ tres .* 1e9; % tres is specified in nanoseconds
stdevs = stdevs ./ tres .* 1e9; % tres is specified in nanoseconds

% Calculate the filtered PSTH
[psthdata_f,zf] = filter(Num,Den,psthdata);
% Compensate for the filter delay
adata = filter(Num,Den,zeros(1,fdelay),zf);
psthdata_f = [psthdata_f(fdelay+1:end) adata];
% Decimate the filtered psth data and contribution, std and stderrs
psthdata_f = psthdata_f(1:dec:end);
contrib = contrib(1:dec:end);
stdevs=stdevs(1:dec:end);
stderrs=stderrs(1:dec:end);
% Calculate the non-filtered, decimated PSTH
if dec > 1
    Numd = ones(1,dec)/dec;
    fdelayd = floor(dec/2);
    [psthdata_nf,zf] = filter(Numd,1,psthdata);
    adata = filter(Numd,1,zeros(1,fdelayd),zf);
    psthdata_nf = [psthdata_nf(fdelayd+1:end) adata];
    psthdata_nf = psthdata_nf(1:dec:end);
else
    psthdata_nf = psthdata;
end

% Calculate the non-contribution-info, decimated PSTH
if dec > 1
    Numd = ones(1,dec)/dec;
    fdelayd = floor(dec/2);
    [psthdata_nco,zf] = filter(Numd,1,psthdata_nco);
    adata = filter(Numd,1,zeros(1,fdelayd),zf);
    psthdata_nco = [psthdata_nco(fdelayd+1:end) adata];
    psthdata_nco = psthdata_nco(1:dec:end);
end

psth.data = psthdata_f;
psth.data_nf = psthdata_nf;
psth.data_nco = psthdata_nco;
psth.stdevs = stdevs;
psth.stderrs = stderrs;
psth.contrib = contrib;
psth.align_offset = floor(align_offset/dec);
psth.time_offset = floor(align_offset/dec)*dec*tres/1e9;
