classdef agc_class < handle
    %agc for a signal of complex points
    %
    
    properties
        agc_val;
        alpha;
        signal_out;
    end
    
    methods
        function Reset(obj)
            obj.agc_val=1;
            obj.alpha=0.1;
            obj.signal_out=[];
        end
        
        function obj = agc_class()
            obj.Reset();
        end
        
        function obj = update(obj,signal_in)
            
%             signal_in(isnan(signal_in))=1;
%             obj.signal_out=signal_in*obj.agc_val;
%             n=numel(signal_in);
%             obj.agc_val=obj.agc_val*(1-obj.alpha)^n+(1-(1-obj.alpha)^n)*(1/mean(abs(signal_in)));
%             return;
            
            obj.signal_out=signal_in;
            %update agc
            for k=1:numel(signal_in)
                if(isnan(signal_in(k)))
                    obj.signal_out(k)=1;
                    continue;
                end
                x=abs(signal_in(k));
                if(obj.agc_val<0.00000001)
                    obj.agc_val=0.00000001;
                end
                if(obj.agc_val>1000000000)
                    obj.agc_val=1000000000;
                end
                obj.agc_val=1/((1/obj.agc_val)*(1-obj.alpha)+obj.alpha*x);
                
                %scale points using agc_val
                obj.signal_out(k)=signal_in(k).*obj.agc_val;
                
            end
        end
    end
    
end

