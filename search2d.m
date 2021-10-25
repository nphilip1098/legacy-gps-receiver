clear all;
close all;
clc;

amount_of_wav_file_to_use_in_seconds=5;
wav_file_offset_in_seconds=0;
c = 2.99792458e8;%speed of light m/s
chips_per_second=1023000;%C/A code chip rate
prn_len=1023;%number of chips till the code starts repeating

prns_per_second=chips_per_second/prn_len;
meters_per_prn=c/prns_per_second;

%read iq samples from wav file and get samplerate
filename='sdrs8.wav';
disp(audioinfo(filename));
[~,Fs] = audioread(filename,[1,1]);%Fs is iq samplerate of wav file
samples = [1+wav_file_offset_in_seconds*Fs,amount_of_wav_file_to_use_in_seconds*Fs+wav_file_offset_in_seconds*Fs];
[y,Fs] = audioread(filename,samples);
y=y(:,1)+1i*y(:,2);

number_of_samples_per_chip=Fs/chips_per_second;
meters_per_sample=c/Fs;

%dont want the 1/2 issue
sympref('HeavisideAtOrigin',1);

sv=17;
prns_per_correlation=10;
prn=2*(cacode(sv,number_of_samples_per_chip)'-0.5);
prn_len_in_samples=numel(prn);
prn_block=repmat(prn,[prns_per_correlation,1]);
prn_block_len_in_samples=numel(prn_block);

signal=y(1:prn_block_len_in_samples);

max_freq_offset_to_try_in_hz=20000;

%calc shift in freq domain for freq offset
hz_per_bin=Fs/prn_block_len_in_samples;
max_freq_shift_to_try=ceil(max_freq_offset_to_try_in_hz/hz_per_bin);

%precomputing for the correlation, basically this is what we are
%doing cor=ifft(fft(signal).*conj(fft(prn))) but the size of a might be
%different so we get frational frequency offsets
A=fft(signal);
cB=conj(fft(prn_block));

%something to save the best correlation and the freq offset of it
maxcorr=0;
maxcorr_freq_shift=0;
maxcorr_chip_shift=0;

%some space for an image
image=zeros(2*max_freq_shift_to_try+1,prn_len_in_samples);

%the 2d search itself
for tmp_freq_shift=-max_freq_shift_to_try:max_freq_shift_to_try
    As=circshift( A, tmp_freq_shift );
    circcorr_xy = ifft(As.*cB);
    image(tmp_freq_shift+max_freq_shift_to_try+1,:)=abs(circcorr_xy(1:prn_len_in_samples));
    [tmp,tmp_chip_shift]=max(abs(circcorr_xy(1:prn_len_in_samples)));% we dont have to go any further than prn_len_in_samples as the prn repeats after that
    if(tmp>maxcorr)
        maxcorr=tmp;
        maxcorr_freq_shift=tmp_freq_shift;
        maxcorr_chip_shift=tmp_chip_shift;
    end
end
surf(image,'linestyle','none','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')

