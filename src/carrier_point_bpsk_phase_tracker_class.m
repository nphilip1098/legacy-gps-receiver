classdef carrier_point_bpsk_phase_tracker_class < handle
    %removes rotation of bpsk points and estimates frequency offset
    %   
    
    properties
        frequency_offset;
        phase_offset;
        
        %for phase error gain
        phase_error_gain;
    end
    
    properties (GetAccess = public , SetAccess = private)
        carrier_phase_error_signal;
    end
    
    properties (GetAccess = private , SetAccess = private)
        point_last;
        fll_b;
        fll_a;
        zi;
    end
    
    properties 
        prn_frequency=1000;%prn freq 1000 prns a second for GPS CA, GAL 250 prns a sec
        
        %for carrier FLL loop filter
        lpf_3db_freq=10;%in hz
        lpf_order=2;%2 should do

    end
    
    methods
        function Reset(obj)
            obj.frequency_offset=0;
            obj.phase_offset=0;
            obj.carrier_phase_error_signal=0;
            obj.point_last=0;
            obj.phase_error_gain=0.18;
            
            %carrier FLL loop filter
            [obj.fll_b,obj.fll_a] = butter(obj.lpf_order,obj.lpf_3db_freq/(obj.prn_frequency/2));
            obj.zi=zeros(obj.lpf_order,1);
            
        end
        
        function set.prn_frequency(obj,value)
            obj.prn_frequency=value;
            %carrier FLL loop filter
            [obj.fll_b,obj.fll_a] = butter(obj.lpf_order,obj.lpf_3db_freq/(obj.prn_frequency/2));
            obj.zi=zeros(obj.lpf_order,1);
        end
        
        function set.lpf_order(obj,value)
            obj.lpf_order=value;
            %carrier FLL loop filter
            [obj.fll_b,obj.fll_a] = butter(obj.lpf_order,obj.lpf_3db_freq/(obj.prn_frequency/2));
            obj.zi=zeros(obj.lpf_order,1);
        end
        
        function set.lpf_3db_freq(obj,value)
            obj.lpf_3db_freq=value;
            %carrier FLL loop filter
            [obj.fll_b,obj.fll_a] = butter(obj.lpf_order,obj.lpf_3db_freq/(obj.prn_frequency/2));
            obj.zi=zeros(obj.lpf_order,1);
        end

        function obj = carrier_point_bpsk_phase_tracker_class()
            obj.Reset();
        end
        
        function points=update(obj,points)
            carrier_freq_error_signal_sum=0;
            for k=1:numel(points)
                point=points(k);
                
                %remove current phase offset
                point_bit=point*exp(-1i*2*pi*obj.phase_offset);
                %calc phase error and update phase_offset
                obj.carrier_phase_error_signal=sign(real(point_bit))*imag(point_bit);
                obj.phase_offset=obj.phase_offset+obj.phase_error_gain*obj.carrier_phase_error_signal;

                %remove nav signal from point before phase offset removed. needed
                %for frequency tracking algo
                if(real(point_bit)<0)
                    point=point*exp(1i*pi);
                end
                
                %calc frequency offset
                carrier_freq_error_signal=imag(point)*real(obj.point_last)-imag(obj.point_last)*real(point);%approx angle between 2 points. delta rad per delta time(1/1000)
                carrier_freq_error_signal=obj.prn_frequency*carrier_freq_error_signal/(2*pi);%convert into Hz
                [carrier_freq_error_signal,obj.zi]=filter(obj.fll_b,obj.fll_a,carrier_freq_error_signal,obj.zi);
                obj.point_last=point;
                carrier_freq_error_signal_sum=carrier_freq_error_signal_sum+carrier_freq_error_signal;
                
                %return the corrected point
                points(k)=point_bit;
            end
            obj.frequency_offset=carrier_freq_error_signal_sum/numel(points);
        end
        
    end
    
end

