function [ output_args ] = ocoher( input_args )
%OCOHER Summary of this function goes here
%   Detailed explanation goes here


end


function [ COH,F,T ] = coher2sigs(sig1,sig2,params)
%COHER2SIGS 2-d coherence basing on matlab spectrogram function
%           for two channels.

% Signal properties
    sig1 = squeeze(sig1);
    sig2 = squeeze(sig2);
    [sigLen nrOfTrials] = size(sig1);
    [sL2 nOT2] = size(sig2);
    assert(sigLen==sL2);
    
    fs = params.fs;
    noverlap = params.no;
    window = params.win;
    nfft = params.nfft;
    k=params.k;
    f=params.f;
    
end