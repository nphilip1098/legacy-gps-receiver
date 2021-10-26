classdef prn_block_nco_class < handle 
    %prn nco that takes a block at a time of M prns
    %   
    
    properties
        phase;         %current phase in chips for this sample (after calling "next", phase will be for the last sample in the vector)
        frequency;     %in chips per ms
        block_len;     %block length in samples
        Fs;            %sample rate in Hz
        sv=1;          %sv number
        gnss;          %GAL or GPS
        use_pilot_prn; %if GAL then you can set this to true and use the pilot prn

    end
    
    properties (GetAccess = public , SetAccess = private)
        prn;
        prn_pilot;
        block;
        
        zero_phase_crossing_index;
        zero_phase_crossing_fractional;
        
        prn_len=1023;
    end
    
    properties (GetAccess = private , SetAccess = private)
        last_phase;
        last_sample_crossing_item;
    end

    
    methods
        
        function Reset(obj)
            obj.use_pilot_prn=false;
            obj.gnss='GAL';
            obj.sv=1;
            obj.Fs=8000000;
            obj.phase=0;
            obj.last_phase=0;
            obj.zero_phase_crossing_index=[];
            obj.last_sample_crossing_item=-inf;
            obj.frequency=2046;
            obj.block_len=(obj.Fs/1000)*obj.prn_len/obj.frequency;%1000 as freq is in ms not s
        end

        function obj = prn_block_nco_class()
            obj.Reset();
        end
        
        function [prn_data,prn_pilot]=gal_calc_prns(obj,sv)
            
            if((sv<0)||(sv>50))
                fprintf("sv out of range\n");
            end
            
            %forget about boc61 as it's too fast

            %create primary codes
            C=GNSScodegen(sv,'E1C');%4092 chips
            B=GNSScodegen(sv,'E1B');%4092 chips
            
            %create subcarrier (only one now)
            boc11=sign(sin(2*pi*(1/2)*[1:2]));%in binary = 10
            
            %expand PRN by subcarrier len
            B2=repelem(B,numel(boc11));
            C2=repelem(C,numel(boc11));
            
            %repeat subcarrier by prn len
            boc11_12=repmat(boc11,[1,numel(B)]);

            %mix PRNs with subcarriers for affective prns
            %we have two prns; take your pick.
            prn_data=(B2.*boc11_12)';%8184 chips
            prn_pilot=(C2.*boc11_12)';%8184 chips

        end

        function set.sv(obj,value)
            obj.sv=value;
            obj.prn=[];
            obj.prn_pilot=[];
            switch(obj.gnss)
                case 'GAL'
                    [obj.prn,obj.prn_pilot]=gal_calc_prns(obj,value);
                    obj.prn_len=8184;
                    obj.frequency=2046;
                case 'GPS'
                    obj.prn=GNSScodegen(value,'L1CA')';
                    obj.prn_len=1023;
                    obj.frequency=1023;
            end
            if(numel(obj.prn)==0)
               error('failed to set prn'); 
            end
        end
        
        function set.gnss(obj,value)
            obj.gnss=value;
            obj.sv=obj.sv;
        end
        
        function prn_block=next(obj)
            
            x=(obj.frequency)*obj.block_len/(obj.Fs/1000);
            a=linspace(0,x,obj.block_len+1)+obj.phase;
            if((obj.use_pilot_prn)&&(strcmp(obj.gnss,'GAL')))
                prn_block=obj.prn_pilot(mod(floor(a(1:obj.block_len)),obj.prn_len)+1);
            else
                prn_block=obj.prn(mod(floor(a(1:obj.block_len)),obj.prn_len)+1);
            end
            obj.block=prn_block;
                
            %for prn crossing
            block_phase=mod((a(1:obj.block_len)),obj.prn_len);
            block_phase_delta=block_phase-[obj.last_phase,block_phase(1:end-1)];
            obj.last_phase=block_phase(end);
            
            obj.zero_phase_crossing_fractional=[];
            obj.zero_phase_crossing_index=[];
            index=find(block_phase_delta<(-obj.prn_len/2));%only interested in +ve crossings of the prn
            for k=1:numel(index)
                delta_samples=index(k)-obj.last_sample_crossing_item;
                obj.last_sample_crossing_item=index(k);
                if(delta_samples>(((obj.Fs/1000)*obj.prn_len/obj.frequency)/2))
%                     fprintf('crossing detected at sample %d of this vector\n',index(k));
                    obj.zero_phase_crossing_index(end+1)=index(k);%integer crossing point
                    obj.zero_phase_crossing_fractional(end+1)=index(k)-block_phase(index(k))/((obj.frequency)/(obj.Fs/1000));%for crossing point with fractional res
                end
            end
            if(numel(index)>0)
                obj.last_sample_crossing_item=index(end)-obj.block_len;
            end
            
            %points to next sample as phase
            %no mod. if you want mod then do that elsewhere
            obj.phase=a(end);

        end
    end
    
end

