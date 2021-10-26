classdef bit_sync_class < handle
    %uses 1khz bpsk points to form 50Hz nav data points and prn index for
    %each 1khz bpsk point
    %

    properties
        bit_sync_matrix_row_size;%number of bits to look over for bit/symbol timing
        
        bit_points;%nav bits (1 is made up of 20 prn points) only valid when prn_number==0.
        prn_numbers;%prn number (0..19 for each prn point)
    end
    
    properties (GetAccess = public , SetAccess = private)
        bit_sync_matrix;
        bit_sync_matrix_row;
        bit_sync_matrix_col;
        bit_sync_matrix_row_start;
        bit_ave;
    end
    
    properties  (Constant)
        bit_sync_matrix_col_size=20;%one bit is 20 points in CA gps
    end
    
    methods
        
        function set.bit_sync_matrix_row_size(obj,value)
            obj.bit_sync_matrix_row_size=value;
            obj.bit_sync_matrix=zeros(obj.bit_sync_matrix_row_size,obj.bit_sync_matrix_col_size);
        end
        
        function Reset(obj)
            obj.bit_sync_matrix_row_size=25;
            obj.bit_sync_matrix_row=1;
            obj.bit_sync_matrix_col=1;
            obj.bit_sync_matrix_row_start=0;
            obj.bit_ave=0;
        end
        
        function obj = bit_sync_class()
            obj.Reset();
        end
        
        function got_bit=bit_sync_update(obj,points)
            got_bit=false;
            obj.bit_points=nan.*zeros(numel(points),1);
            obj.prn_numbers=zeros(numel(points),1);
            for n=1:numel(points)
                point=points(n);
                
                %put point into bit sync matrix
                obj.bit_sync_matrix_col=mod(obj.bit_sync_matrix_col,obj.bit_sync_matrix_col_size)+1;
                if(obj.bit_sync_matrix_col==1)
                    obj.bit_sync_matrix_row=mod(obj.bit_sync_matrix_row,obj.bit_sync_matrix_row_size)+1;
                end
                obj.bit_sync_matrix(obj.bit_sync_matrix_row,obj.bit_sync_matrix_col)=heaviside(real(point));
                
                %find col with most bit transitions
                val_max=0;
                val_max_k_next=1;
                for k=1:obj.bit_sync_matrix_col_size
                    k_next=mod(k,obj.bit_sync_matrix_col_size)+1;
                    if(k_next==1)
                        val=sum(bitxor(obj.bit_sync_matrix(1:end-1,k),obj.bit_sync_matrix(2:end,k_next)));
                    else
                        val=sum(bitxor(obj.bit_sync_matrix(:,k),obj.bit_sync_matrix(:,k_next)));
                    end
                    if(val>val_max)
                        val_max=val;
                        val_max_k_next=k_next;
                    end
                end
                if(obj.bit_sync_matrix_row_start~=val_max_k_next)
                    fprintf('bit sync change was %d now %d\n',obj.bit_sync_matrix_row_start,val_max_k_next);
                end
                obj.bit_sync_matrix_row_start=val_max_k_next;
                
                %when point of start of bit arrives
                if(obj.bit_sync_matrix_row_start==obj.bit_sync_matrix_col)
                    obj.bit_ave=obj.bit_ave./obj.bit_sync_matrix_col_size;
                    
                    %add bit to return. we probably will only get one bit
                    %at most per call but anyway.
                    got_bit=true;
                    obj.bit_points(n)=obj.bit_ave;
                    
                    obj.bit_ave=0;
                end
                obj.bit_ave=obj.bit_ave+point;
                                
                %add return of what prn number this is 0..19 
                obj.prn_numbers(n)=mod(obj.bit_sync_matrix_col-obj.bit_sync_matrix_row_start,obj.bit_sync_matrix_col_size);
                
            end
        end
    end

end