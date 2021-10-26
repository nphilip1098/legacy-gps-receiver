classdef search_2d_for_prn_class < handle
    %search for prn in chip phase and cariier frequency
    %   
    
    properties
        prn_block;%consists of many prns
        prn_len;%length of 1 prn in the prn_block
        max_freq_offset_to_try_in_hz;
        Fs;
    end
    
    properties (GetAccess = public , SetAccess = private)
        frequency_offset_in_hz;
        chip_offset_in_samples;
        correlation;
    end
    
    properties (GetAccess = private , SetAccess = private)
        cB;
    end
    
    methods
        
        function Reset(obj)
            obj.max_freq_offset_to_try_in_hz=10000;
            obj.Fs=8000000;
            obj.prn_len=8000;
            obj.prn_block=zeros(40000,1);
            obj.frequency_offset_in_hz=0;
            obj.chip_offset_in_samples=0;
            obj.correlation=1;
        end

        function obj = search_2d_for_prn_class()
            obj.Reset();
        end
        
        function set.prn_block(obj,value)
            obj.prn_block=value;
            %precomputing for the correlation, basically this is what we are
            %doing cor=ifft(fft(signal).*conj(fft(prn))) but the size of a might be
            %different so we get frational frequency offsets
            obj.cB=conj(fft(obj.prn_block));
        end
        
        function obj = search(obj,signal)
            
            prn_block_len_in_samples=numel(obj.prn_block);
            
            assert(prn_block_len_in_samples==numel(signal),'signal and prn block should be same lengths');
            assert(mod(prn_block_len_in_samples,obj.prn_len)==0,'number of prns in prn_block needs to be a whole number');
            
            %calc shift in freq domain for freq offset
            hz_per_bin=obj.Fs/prn_block_len_in_samples;
            max_freq_shift_to_try=ceil(obj.max_freq_offset_to_try_in_hz/hz_per_bin);
            assert(max_freq_shift_to_try>=0,'max_freq_shift_to_try should not be negative');

            %precomputing for the correlation, basically this is what we are
            %doing cor=ifft(fft(signal).*conj(fft(prn))) but the size of a might be
            %different so we get frational frequency offsets
            A=fft(signal);

            %something to save the best correlation and the freq offset of it
            maxcorr=0;
            maxcorr_freq_shift=0;
            maxcorr_chip_shift=0;
% % % % freqimage=[];            
% % % % % image=zeros(2*max_freq_shift_to_try+1,obj.prn_len);

            %the 2d search itself
            for tmp_freq_shift=-max_freq_shift_to_try:max_freq_shift_to_try
                As=circshift( A, tmp_freq_shift );
                circcorr_xy = ifft(As.*obj.cB);
% % % % % image(tmp_freq_shift+max_freq_shift_to_try+1,:)=abs(circcorr_xy(1:obj.prn_len));
                [tmp,tmp_chip_shift]=max(abs(circcorr_xy(1:obj.prn_len)));% we dont have to go any further than prn_len_in_samples as the prn repeats after that
% % % % freqimage(end+1)=tmp;
                if(tmp>maxcorr)
                    maxcorr=tmp;
                    maxcorr_freq_shift=tmp_freq_shift;
                    maxcorr_chip_shift=tmp_chip_shift;
                end
            end
            
% % % % % %             As=circshift( A, maxcorr_freq_shift );
% % % % % %             circcorr_xy = ifft(As.*obj.cB);
% % % % % %             plot(abs(circcorr_xy));
% % % % % %             jbjhb


% plot((-max_freq_shift_to_try:max_freq_shift_to_try)*hz_per_bin-maxcorr_freq_shift*hz_per_bin,freqimage);
% xlabel('frequency offset from max max frequency (Hz)');
% % % % plot(round((-max_freq_shift_to_try:max_freq_shift_to_try)*hz_per_bin),freqimage);
% % % % xlabel('frequency offset (Hz)');
% % % % figure;
% % % % 
% % % % hz_per_bin
% asdasda
            
% % % % % surf(image,'linestyle','none','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')
% % % % % ccghgch
            
            %calc best freq, chip offsets, and correlation
            obj.frequency_offset_in_hz=maxcorr_freq_shift*hz_per_bin;
            obj.chip_offset_in_samples=maxcorr_chip_shift;
            obj.correlation=maxcorr/prn_block_len_in_samples;
            
        end
    end
    
end

