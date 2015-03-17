function [ coherenceMaps,freqVec,timeVec ] = ocoher(signal,params)
%OCOHER 2-d coherence basing on matlab spectrogram function.
% signal is a cell with experimental conditions and each cell
% is a matrix of size: (nrOfChannels signalLength nrOfTrials)

    nrOfConditions = length(signal);
    [nrOfChannels sigLen nrOfTr] = size(signal{1,1});

        % Spectrogram parameters
    if ~isfield(params,'fs')
        params.fs = 512;
    end
    if ~isfield(params,'win')
        params.win=128;
    end
    if ~isfield(params,'no')
        params.no = params.win-2;
    end
    if ~isfield(params,'nfft')
        params.nfft = 256;
    end

    %spectrogram size
    k = fix((sigLen-params.no)/(params.win-params.no));
    if mod(params.nfft,2)==0
        f = (params.nfft/2+1);
    else
        f = (params.nfft+1)/2;
    end
    params.k = k;
    params.f = f;
    coherenceMaps = cell(1,nrOfConditions);

    for conditionNr = 1:nrOfConditions
        coherenceMaps{1,conditionNr} = zeros(f,k,nrOfChannels,nrOfChannels);
    end

    
    for conditionNr = 1:nrOfConditions
        for chan1 = 1:nrOfChannels        
            for chan2 = chan1 : nrOfChannels
                disp([' conditionNumber: ',num2str(conditionNr),', channels: ',num2str(chan1),' -> ',num2str(chan2)]);
                [COH,F,T] = coher2sigs(signal{conditionNr}(chan1,:,:),signal{conditionNr}(chan2,:,:),params);
                coherenceMaps{conditionNr}(:,:,chan1,chan2) = COH;
                coherenceMaps{conditionNr}(:,:,chan2,chan1) = -1*COH;
            end
        end
    end
    freqVec = F; timeVec = T;
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
    
    %Cross spectrum calculating:
    crossspect = zeros(f,k);
    cfs_S1 = zeros(f,k);
    cfs_S2 = zeros(f,k);
    for trNr = 1 : nrOfTrials
        [S1,F,T]   = spectrogram(sig1(:,trNr),window,noverlap,nfft,fs);
        [S2,F,T]   = spectrogram(sig2(:,trNr),window,noverlap,nfft,fs);
        cfs        = S1.*conj(S2);
        crossspect = crossspect + cfs;
        cfs_S1     = cfs_S1 + abs(S1);
        cfs_S2     = cfs_S2 + abs(S2);        
    end
    cfs_S1     = cfs_S1/nrOfTrials;
    cfs_S2     = cfs_S2/nrOfTrials;
    crossspect = crossspect/nrOfTrials;
    
    COH = crossspect./(cfs_S1.*cfs_S2);

end