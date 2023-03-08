function [meas] = FilterBoomSignals(meas,dt,fc,fo,filter_type)
%FilterBoomSignals low-pass filters the boom air data signals

% Retrieve boom signals from struct
alpha = meas(11,:);
beta = meas(12,:);

% Apply the low-pass filter
alpha_filt = LPfilter(alpha, fc, fo, dt, filter_type);
beta_filt = LPfilter(beta, fc, fo, dt, filter_type);

% Store the filtered data in 
meas(11,:) = alpha_filt;
meas(12,:) = beta_filt;
end