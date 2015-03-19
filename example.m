%signal preparing:
fs = 1024;
nr = 20; % number of trials
timeLen = 1;
sig(1,:,:) =  zeros(timeLen*fs,nr);
sig(2,:,:) =  zeros(timeLen*fs,nr);
t   = linspace(0,timeLen,timeLen*fs);
fq = 34;
lb = -4; ub =4;
[psi,x] = mexihat(lb,ub,timeLen*fs);
delay = fix(0.5*timeLen*fs);
for i=1:nr
    sig(1,:,i) = sin(2*pi*fq*t)+0.6*randn(1,timeLen*fs);
    sig(1,1:delay,i) = 0.3*randn(1,delay);
    sig(2,:,i) = 5*psi.*sin(2*pi*fq*t)+0.6*randn(1,timeLen*fs);
end

%coherence calculating:

params.fs = fs;
params.nfft = 512;
params.no = 255;
params.win = 256;

[C,F,T] = ocoher({sig},params);
imagesc(T,F,abs(C{1}(:,:,1,2)));set(gca,'Ydir','normal');