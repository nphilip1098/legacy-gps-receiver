clear all;
close all;
clc;
format long g;


%wav filename to use
%only tested with 8Ms/s samplerate

filename='baseband_1574740000Hz_09-17-12_26-10-2021.wav';

%amount and where to start reading the wav file
amount_of_wav_file_to_use_in_seconds=5;
wav_file_offset_in_seconds=0;

%this is set in srduno. currently only needed for carrier aiding
tuner_frequency_offset_for_carrier_aiding=0;

%carrier aiding is less noisy for the prn tracking
use_carrier_aiding=true;

%either use GAL or GPS
gnss_system='GPS';

%search for visible satelites
search_for_detectable=true;

%stop if doing just a search
stop_once_search_for_detectable_finished=true;
  
%svs to use if search_for_detectable is false
% if(strcmp(gnss_system,'GPS'))
%         svs=[ 12 6 17 19 3 14 28 24 2];%gps for demo.wav
% else
%         svs=[ 25 12 2 11 14 24];%gal for demo.wav
% end    

%index of svs to use
sv_index=1;

%if GAL then setting this will cause tracking to be done on the pilot
gal_use_pilot=true;

%only usfull if amount_of_wav_file_to_use_in_seconds is small say 2 seconds
%or so
gal_plot_mod25_demod=true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% things should be setup for searching and tracking from here %
% may need changing if sample rate isn't 8Ms/s                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 2.99792458e8;%speed of light m/s
chips_per_second=1023000;%chip rate
prn_length=1023;%number of chips till the code starts repeating
if(strcmp(gnss_system,'GAL'))
    chips_per_second=1023000*2;%effective boc11 code chip rate
    prn_length=1023*4*2;%number of boc11 chips till the code starts repeating
end
prns_per_second=chips_per_second/prn_length;% prns/s
meters_per_prn=c/prns_per_second;%

disp(audioinfo(filename))
[~,Fs] = audioread(filename,[1,1]);%Fs is iq samplerate of wav file
samples = [1+wav_file_offset_in_seconds*Fs,amount_of_wav_file_to_use_in_seconds*Fs+wav_file_offset_in_seconds*Fs];
disp(samples);
[y,Fs] = audioread(filename,samples);
y=y(:,1)+1i*y(:,2);
% disp(y);

number_of_samples_per_chip=Fs/chips_per_second;
meters_per_sample=c/Fs;

if(search_for_detectable)
    
    %number of wanted prns to correlate over for searching
    prns_per_correlation=10;
    if(strcmp(gnss_system,'GAL'))
        prns_per_correlation=1;
    end
    
    samples_per_prn=Fs/prns_per_second;
    samples_per_correlation=samples_per_prn*prns_per_correlation;
    
    % prn nco
    prn_nco=prn_block_nco_class();
    prn_nco.gnss=gnss_system;
    prn_nco.block_len=samples_per_correlation;
    prn_nco.Fs=Fs;
    prn_nco.sv=1;
    if(strcmp(gnss_system,'GAL'))
        prn_nco.use_pilot_prn=gal_use_pilot;
    end
    %
    
    % 2d search
    %GAL has 2 paths from the prns to the real output adn we dont know what the
    %secondry prn number or the INAV data is and these mix with the PRNs
    search_2d=search_2d_for_prn_class();
    search_2d.Fs=Fs;
    search_2d.max_freq_offset_to_try_in_hz=20000;
    search_2d.prn_block=prn_nco.next();
    search_2d.prn_len=samples_per_prn;
    %
    
    max_sv_id=32;
    if(strcmp(gnss_system,'GAL'))
        max_sv_id=36;
    end
    svs_corr=zeros(max_sv_id,1);
    for sv=1:max_sv_id
        found_sv=false;
        
        %set sv number to try
        prn_nco.sv=sv;
        search_2d.prn_block=prn_nco.next();
        
        %look through some blocks and take ave.
        %we know we can get a big
        %reduction in corr so this will help
        %to avoid that issue.
        abs_corr_sum=0;
        n_tries=6;
        if(strcmp(gnss_system,'GAL'))
            n_tries=2;
        end
        for m=0:(n_tries-1)
            %get a block of signal
            a_org=y(samples_per_correlation*m+1:samples_per_correlation*(m+1));
            %search
            search_2d.search(a_org);
            %add to sum
            abs_corr_sum=abs_corr_sum+search_2d.correlation;
        end
        
        %print and save corr to array
        fprintf("sv=%d corr=%f\n",sv,abs_corr_sum);
        svs_corr(sv)=abs_corr_sum;
    end
    
    %use built in function to find the outliers (detected signals)
    svs=find(isoutlier(svs_corr));
    %think this is the same
    %as octave doesn't have this function
    ac=-1/(sqrt(2)*erfcinv(3/2));
    svs=find(abs(svs_corr-median(svs_corr))>3*ac*median(abs(svs_corr-median(svs_corr))));
    
    %sort in order
    [~,idx]=sort(svs_corr(svs),'descend');
    svs=svs(idx);
    
    fprintf('prns are the following. copy this into the script\n');
    fprintf('svs=[ ');fprintf("%d ",svs');fprintf('];\n');
    
    if(stop_once_search_for_detectable_finished)
      return;
    end

end

sv=svs(sv_index);
results_for_ls_solving(sv_index).sv=sv;
prns_per_correlation=15;
if(strcmp(gnss_system,'GAL'))
    prns_per_correlation=2;
end

samples_per_prn=Fs/prns_per_second;
samples_per_correlation=samples_per_prn*prns_per_correlation;

%prn nco
prn_nco=prn_block_nco_class();
prn_nco.gnss=gnss_system;
prn_nco.block_len=samples_per_correlation;
prn_nco.Fs=Fs;
prn_nco.sv=sv;
if(strcmp(gnss_system,'GAL'))
    prn_nco.use_pilot_prn=gal_use_pilot;
end
%

%gains for EL prn_nco freq and phase update
chip_tracking_error_phase_gain=0.25/2;
chip_tracking_error_freq_gain=(1/2000)*2;

%integration and dump
integrate_and_dump=integrate_and_dump_class();
%

%agc
agc=agc_class();
%

%phase tracker
carrier_tracker=carrier_point_bpsk_phase_tracker_class();
carrier_tracker.phase_error_gain=0.18/2;
if(strcmp(gnss_system,'GAL'))
    carrier_tracker.phase_error_gain=0.18/2;
end
%

%2d search
search_2d=search_2d_for_prn_class();
search_2d.Fs=Fs;
search_2d.max_freq_offset_to_try_in_hz=20000;
search_2d.prn_block=prn_nco.next();
search_2d.prn_len=samples_per_prn;
%

%for GPS. sub sampling not needed as GPS has 500Hz bandwidth
sub_sample_factor=1;

%for GAL's narrowband of 125Hz rather than GPS's 500Hz
if(strcmp(gnss_system,'GAL'))
    carrier_tracker_fast_points=carrier_point_bpsk_phase_tracker_class();
    carrier_tracker_fast_points.phase_error_gain=0.18/4;
    sub_sample_factor=4;
    carrier_tracker_fast_points.prn_frequency=250*sub_sample_factor;
    carrier_tracker_fast_points.lpf_3db_freq=125;
end

%logs
prn_point_log=[];
agc_log=[];
freq_offset_log=[];
chip_tracking_error_log=[];
carrier_phase_error_signal_log=[];
prn_freq_offset_log=[];
time_log=[];

%for carrier freq that has no class for it
frequency_offset_gain=0.25*2;%might be a be too big but locks fast
if(strcmp(gnss_system,'GAL'))
    frequency_offset_gain=0.35;
end
freq_offset=0;
carrier_freq_phase_remainder=0;
carrier_freq_phase_count=0;

%using the aparent 1.023*2 chip rate for the sinboc(1,1)
carrier_aiding_factor=1;
if(strcmp(gnss_system,'GAL'))
    carrier_aiding_factor=2;
end

%for switching btween searching and tracking states
locked=false;

%display settings
fprintf('--------------------\n');
fprintf('GNSS system %s\n',gnss_system)
fprintf('sv=%d\n',sv);
fprintf('prns_per_correlation=%d\n',prns_per_correlation);
fprintf('--for carrier tracking--\n');
fprintf('carrier_tracker.phase_error_gain=%f\n',carrier_tracker.phase_error_gain);
fprintf('frequency_offset_gain=%f\n',frequency_offset_gain);
fprintf('--for prn tracking--\n');
fprintf('chip_tracking_error_phase_gain=%f\n',chip_tracking_error_phase_gain);
fprintf('chip_tracking_error_freq_gain=%f\n',chip_tracking_error_freq_gain);
fprintf('--------------------\n');

%start timer
tic

%run through blocks of signal
for m=1:floor((numel(y)-1)/samples_per_correlation-1)
    
    %approx start time in secs
    block_start_sample=samples_per_correlation*m+1;
    block_approx_start_time_wrt_start_of_signal=block_start_sample/Fs;
    
    %get a block of prn from the prn nco
    prn_nco.next();
    %get a block of signal
    a_org=y(samples_per_correlation*m+1:samples_per_correlation*(m+1));

    %if not locked
    if(~locked)

        %search for signal changing the prn phase given a nominal prn freq
        %chips/s chip rate and also the carrier freq
        search_2d.search(a_org);
        %use the search result to adjust the carrier freq,
        %the prn nco phase, and also the agc gain.
        freq_offset=search_2d.frequency_offset_in_hz;
        prn_nco.phase=-(search_2d.chip_offset_in_samples-1)/number_of_samples_per_chip;
        prn_nco.next();
        agc.agc_val=(1/search_2d.correlation);
        carrier_freq_phase_remainder=0;
        
        carrier_tracker.phase_offset=0;
        carrier_tracker.frequency_offset=0;
        
        carrier_tracker_fast_points.phase_offset=0;
        carrier_tracker_fast_points.frequency_offset=0;
        
        carrier_freq_phase_count=0;
        
        %print the search result
        fprintf('--2d search estimation--\n');
        fprintf('frequency_offset_in_hz=%f\n',freq_offset);
        fprintf('prn_nco_phase=%f\n',prn_nco.phase);
        fprintf('agc_val=%f\n',agc.agc_val);
        
        %say we are locked
        locked=true;
        
    end
    
    %carrier_freq_phase,carrier_freq_phase_remainder and
    %carrier_freq_phase_count are for last sample in a_org
    %after the floowoing is done.
    a_slow_spinning=a_org.*exp( ...
        -2*pi*1i*([1:samples_per_correlation].*(freq_offset/Fs)+carrier_freq_phase_remainder) ...
        )';
    carrier_freq_phase=samples_per_correlation*(freq_offset/Fs)+carrier_freq_phase_remainder;
    carrier_freq_phase_remainder=mod(carrier_freq_phase,1);
    carrier_freq_phase_count=carrier_freq_phase_count+round(carrier_freq_phase-carrier_freq_phase_remainder);
    
    %correlate slow spinning signal with prn block and integrate and dump
    %over each prn not crossing prn boundries
    integrate_and_dump.next(a_slow_spinning.*prn_nco.block,prn_nco.zero_phase_crossing_index);
    %agc for normalized dumped points
    agc.update(integrate_and_dump.dumps./integrate_and_dump.dumps_count);
    
    %carrier phase tracking
    agc.signal_out=carrier_tracker.update(agc.signal_out);
    
    %narrow band GAL
    %this will sub sample a prn by sub_sample_factor times so that
    %carrier_tracker_fast_points can see if the point is moving.
    %this is handy for GAL where prns are 4ms giving a locking bandwidth of
    %+-62.5Hz. if sub_sample_factor=8 then becomes +-500Hz witch the FFT
    %will find easy.
    if(sub_sample_factor>1)
        fast_points=agc.agc_val*sum(reshape(a_slow_spinning.*prn_nco.block,[samples_per_correlation/(sub_sample_factor*prns_per_correlation),(sub_sample_factor*prns_per_correlation)]))/(samples_per_correlation/(sub_sample_factor*prns_per_correlation));
        fast_points=carrier_tracker_fast_points.update(fast_points);
        carrier_tracker.frequency_offset=carrier_tracker_fast_points.frequency_offset;
    end
    
    %update carrier frequency (not phase!)
    %unlock if it goes out of bounds
    freq_offset=freq_offset-frequency_offset_gain*carrier_tracker.frequency_offset;

    if(abs(freq_offset)>search_2d.max_freq_offset_to_try_in_hz)
        locked=false;
        fprintf('unlocked at %.3fs due to carrier offset too big\n',block_approx_start_time_wrt_start_of_signal);
    end
    
    %unlock if prn corr goes too low
    for yy=1:numel(agc.signal_out)
        if(abs(agc.signal_out(yy))<0.125)
            locked=false;
            fprintf('unlocked at %.3fs\n',block_approx_start_time_wrt_start_of_signal);
        end
    end
    
    %calc EL chip error
    matched_prn_early=circshift( prn_nco.block, -1 );
    matched_prn_late=circshift( prn_nco.block, +1 );
    chip_tracking_early=(sum(a_slow_spinning.*matched_prn_early)*(agc.agc_val/samples_per_correlation));
    chip_tracking_late=(sum(a_slow_spinning.*matched_prn_late)*(agc.agc_val/samples_per_correlation));
    chip_tracking_error=(abs(chip_tracking_early)-abs(chip_tracking_late));

    %use EL chip error to adjust chip phase. frequency uses either
    %carrier aiding or EL chip error
    prn_nco.phase=prn_nco.phase+chip_tracking_error*chip_tracking_error_phase_gain/2;
    %if carrier aiding
    if(use_carrier_aiding)
        prn_nco.frequency=((tuner_frequency_offset_for_carrier_aiding-freq_offset)/(1540*1000)+1023)*carrier_aiding_factor;
    else
        prn_nco.frequency=prn_nco.frequency+chip_tracking_error*chip_tracking_error_freq_gain;
    end

    carrier_phase_error_signal_log(end+1)=carrier_tracker.carrier_phase_error_signal;
    agc_log(end+1)=agc.agc_val;
    chip_tracking_error_log(end+1)=chip_tracking_error;
    freq_offset_log(end+1)=freq_offset;
    prn_freq_offset_log(end+1)=prn_nco.frequency;
    prn_point_log=cat(1,prn_point_log,agc.signal_out');
    time_log(end+1)=samples_per_correlation*(m+1)/Fs;
    
end
toc

plot(prn_point_log(25:end),'.')%you might want to remove the first bit of the signal just to make it clearer
xlim([-2,2]);
ylim([-2,2]);
grid on;
title({[ 'PRN ' num2str(sv) ' correlation points']});

%this doesn't work in octave. seems the small ydiff is the problem
figure;
plot(time_log,prn_freq_offset_log)
grid on;
title({[ 'PRN ' num2str(sv) ' chip frequency']});
ylabel('Frequency (kHz)');
xlabel('Time (s)');

figure;
plot(time_log,freq_offset_log)
title({[ 'PRN ' num2str(sv) ' carrier freq offset']});
ylabel('Frequency (Hz)');
xlabel('Time (s)');

figure;
plot(time_log,agc_log);
title({[ 'PRN ' num2str(sv) ' agc val']});
xlabel('Time (s)');

figure;
plot(time_log,chip_tracking_error_log);
title({[ 'PRN ' num2str(sv) ' chip tracking error']});
ylabel('Error (???)');
xlabel('Time (s)');

figure;
plot(time_log,carrier_phase_error_signal_log);
title({[ 'PRN ' num2str(sv) ' carrier tracking error']});
ylabel('Error (???)');
xlabel('Time (s)');

%GAL only stuff
if(strcmp(gnss_system,'GAL'))
    
    
    mod_val=25;
    
    rx_data=heaviside(real(prn_point_log(1:end)),1);
    
    if(gal_use_pilot)
        
        %find pilot signal and align and invert if needed
        number_of_start_bits_to_search=75;
        number_of_tests=5;
        assert(((number_of_start_bits_to_search-1)+(25*number_of_tests))<=numel(rx_data),'not enough demodulated bits\n');
        secondary_prn=[0 0 1 1 1 0 0 0 0 0 0 0 1 0 1 0 1 1 0 1 1 0 0 1 0 ]';%CS25
        secondary_prn=repmat(secondary_prn,[1,number_of_tests]);
        val_min=inf;
        val_max=-inf;
        min_loc=-1;
        max_loc=-1;        
        for k=1:number_of_start_bits_to_search
            tmp=reshape(rx_data(k:((k-1)+25*number_of_tests)),[25,number_of_tests]);
            val=sum(sum(bitxor(secondary_prn,tmp))/(25*number_of_tests));
            if(val<val_min)
                val_min=val;
                min_loc=k;
            end
            if(val>val_max)
                val_max=val;
                max_loc=k;
            end
        end
        if((1-val_max)<val_min)
            fprintf('inverted signal\n');
            rx_data=1-rx_data;
            loc=max_loc;
        else
            loc=min_loc;
        end
        loc=mod(loc-1,25)+1;
        rx_data=rx_data(loc:end);
        
    end
    
    fprintf('some data in 25 bit blocks\n');
    for n=0:floor(numel(rx_data)/25-1)
        strvar="";
        for k=1:25
            strvar=strcat(strvar,int2str(rx_data(k+n*25)));
        end
        fprintf("%s\n",strvar);
    end
    
    if(gal_plot_mod25_demod)
        figure;
        rx_data=rx_data(1:mod_val*floor(numel(rx_data)/mod_val));
        rx_data=rx_data+1.1*floor([0:numel(rx_data)-1]'/mod_val);
        rx_data = reshape(rx_data,[mod_val,numel(rx_data)/mod_val]);
        plot([0:mod_val-1],rx_data);
        xlim([0 mod_val-1]);
        ylim([-0.1 size(rx_data,2)*1.1]);
        xlabel(['PRN mod ' num2str(mod_val)]);
        ylabel('Demod trace');
        title({['Demod trace vs PRN ' num2str(sv) ' mod ' num2str(mod_val)]});
    end

end

