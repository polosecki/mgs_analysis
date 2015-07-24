function [psthfilter] = design_psth_filter(psthfilter)
% Compute filter coefficients based on the psthfilter parameters
% struct.

% Design the filter
if strcmp(psthfilter.type,'Gauss')
    fsms = round(psthfilter.fs/1000);
    sigmas = psthfilter.sigma_ms/1000;
    x = [-psthfilter.length_ms*fsms:-1 0 1:psthfilter.length_ms*fsms]/psthfilter.fs;
    Num = 1 ./ (sqrt(2*pi)*sigmas) .* exp(-x.^2 / (2*sigmas^2));
    Num = Num/psthfilter.fs;
    psthfilter.Num = Num;
    psthfilter.Den = 1;
    psthfilter.delay = round((length(Num)-1)/2);
elseif strcmp(psthfilter.type,'Dolph-Chebyshev')
    cd dolph_chebyshev
    Num = dolph_chebyshev(psthfilter.fs,psthfilter.fc,psthfilter.as);
    cd ..
    psthfilter.Num = Num;
    psthfilter.Den = 1;
    psthfilter.delay = round((length(Num)-1)/2);
elseif strcmp(psthfilter.type,'Butterworth')
    [psthfilter.Num,psthfilter.Den] = butter(psthfilter.order,psthfilter.fc*2/psthfilter.fs);
    psthfilter.delay = 0;
else
    error('You specified an unknown PSTH filter type');
end
