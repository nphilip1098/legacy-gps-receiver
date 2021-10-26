
clear all;
clc;

amount_of_wav_file_to_use_in_seconds=1;
wav_file_offset_in_seconds=0;
c = 2.99792458e8;%speed of light m/s
chips_per_second=1023000;%C/A code chip rate
prn_len=1023;%number of chips till the code starts repeating

prns_per_second=chips_per_second/prn_len;
meters_per_prn=c/prns_per_second;

%read iq samples from wav file and get samplerate
filename='demo.wav';
[~,Fs] = audioread(filename,[1,1]);%Fs is iq samplerate of wav file
samples = [1+wav_file_offset_in_seconds*Fs,amount_of_wav_file_to_use_in_seconds*Fs+wav_file_offset_in_seconds*Fs];
[y,Fs] = audioread(filename,samples);
y=y(:,1)+1i*y(:,2);

number_of_samples_per_chip=Fs/chips_per_second;
meters_per_sample=c/Fs;

%dont want the 1/2 issue
sympref('HeavisideAtOrigin',1);

% find what svs we can see
%number of wanted prns to correlate over

prns_per_correlation=20;

samples_per_prn=Fs/prns_per_second;
samples_per_correlation=samples_per_prn*prns_per_correlation;

% prn nco
prn_nco=prn_block_nco_class();
prn_nco.block_len=samples_per_correlation;
prn_nco.Fs=Fs;
prn_nco.sv=1;
prn_nco.phase=0;
prn_nco.frequency=1023;
%

% 2d search
search_2d=search_2d_for_prn_class();
search_2d.Fs=Fs;
search_2d.max_freq_offset_to_try_in_hz=20000;
search_2d.prn_block=prn_nco.next();
search_2d.prn_len=samples_per_prn;
%

svs_corr=zeros(32,1);
for sv=1:32
    found_sv=false;
    
    %set sv number to try
    prn_nco.sv=sv;
    search_2d.prn_block=prn_nco.next();
    
    %look through 4 blocks and take max.
    %we know we can get a big
    %reduction in corr so this will help
    %to avoid that issue.
    abs_corr_sum=0;
    for m=0:3
        %get a block of signal
        a_org=y(samples_per_correlation*m+1:samples_per_correlation*(m+1));
        
        %search for signal changing the prn phase given a nominal 1023
        %chips/s chip rate and also the carrier freq
        search_2d.search(a_org);
        
        abs_corr_sum=abs_corr_sum+search_2d.correlation;
    end
    
    %print and save corr to array
    fprintf("sv=%d corr=%f\n",sv,abs_corr_sum);
    svs_corr(sv)=abs_corr_sum;
end

%sort corr in decending order
[svs_corr_sorted,svs_sorted]=sort(svs_corr,'descend');
%we can only see 50% of sv max on earth so use the weakest svs as noise ref
%and remove this from our result. anything above we assume to be a signal
svs_corr_sorted=svs_corr_sorted-max(svs_corr_sorted(17:end))-3*std(svs_corr_sorted(17:end));
svs_corr_sorted=svs_corr_sorted./svs_corr_sorted(1);
%display the results
plot(svs_sorted,svs_corr_sorted,'o')
ylim([0 svs_corr_sorted(1)+3*std(svs_corr_sorted(17:end))./svs_corr_sorted(1)]);
xlim([0 33]);
%make the a list of detectable svs
svs=svs_sorted(svs_corr_sorted>0);