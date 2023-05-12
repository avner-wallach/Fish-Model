function [ Zout ] = Trans2LFP(Zin,params)
% TRANS2LFP convert model's 'transdermal potential' into 'LFP' using
% non-linear function
th=params.th;
sig=params.sig;

Zout=tanh((Zin-th)/sig);


end

