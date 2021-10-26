classdef integrate_and_dump_class < handle
    %intgrate and dump a signal
    %
    
    properties
        dumps;
        dumps_count;
    end
    
    properties (GetAccess = private , SetAccess = private)
        running_integration_sum;
        running_integration_count;
    end
    
    methods
        
        function Reset(obj)
            obj.running_integration_sum=0;
            obj.running_integration_count=0;
            obj.dumps=[];
            obj.dumps_count=[];
        end
        
        function obj = integrate_and_dump_class()
            obj.Reset();
        end
        
        function obj = next(obj,signal,dump_index)
            obj.dumps=[];
            obj.dumps_count=[];
            last_sample=1;
            for k=1:numel(dump_index)
                start_sample=dump_index(k);
                this_count=start_sample-1-last_sample+1;
                this_sum=sum(signal(last_sample:start_sample-1));
                obj.running_integration_count=obj.running_integration_count+this_count;
                obj.running_integration_sum=obj.running_integration_sum+this_sum;
                
                %dump the sums for the user
                obj.dumps(end+1)=obj.running_integration_sum;
                obj.dumps_count(end+1)=obj.running_integration_count;
                
%                 if(last_sample==1)
%                     fprintf('first dump: sum=%f count=%d\n',obj.running_integration_sum,obj.running_integration_count);
%                 else
%                     fprintf('middle dump: sum=%f count=%d\n',obj.running_integration_sum,obj.running_integration_count);
%                 end
                
                last_sample=start_sample;
                obj.running_integration_sum=0;
                obj.running_integration_count=0;
            end
            
            %save partial sums for later
            this_count=numel(signal)-last_sample+1;
            this_sum=sum(signal(last_sample:end));
            obj.running_integration_count=obj.running_integration_count+this_count;
            obj.running_integration_sum=obj.running_integration_sum+this_sum;
            
%             if(obj.running_integration_count~=0)
%                 fprintf('end store: sum=%f count=%d\n',obj.running_integration_sum,obj.running_integration_count);
%             end
            
        end
        
        
    end
    
end

