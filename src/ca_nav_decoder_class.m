classdef ca_nav_decoder_class < handle
    %decodes c/a nav data
    %

    properties
        tow;
        subframe;
        subframe_bad_word_count;%if not 10 then at least TML and HOW are valid
    end
    
    properties (GetAccess = private , SetAccess = private)
        bit_sub_frame_buffer;
    end
    
    properties  (Constant)
        
    end
    
    methods
        
        function Reset(obj)
            obj.bit_sub_frame_buffer=zeros(1,2+10*30);%2 bits for D29 and D30 of last word and 10 words of length 30 for the data
            obj.subframe=zeros(1,10);%to return words
            obj.subframe_bad_word_count=10;
            obj.tow=nan;
        end
        
        function obj = ca_nav_decoder_class()
            obj.Reset();
        end
        
        %50Hz points one at a time please
        function update(obj,bit_point)
            
            obj.subframe_bad_word_count=10;
            obj.tow=nan;
    
            %shuffle things over and put in the new bit
            obj.bit_sub_frame_buffer=circshift(obj.bit_sub_frame_buffer,-1);
            obj.bit_sub_frame_buffer(end)=heaviside(real(bit_point));
            
            %see if we have a good TLM and HOW
            [ok,TLM_word]=paritycheck(obj.bit_sub_frame_buffer(1+30*(1-1):32+30*(1-1)));
            TLM_preamble=bitand(bitshift(TLM_word,-16),255);
            if((ok)&&(TLM_preamble==139))%bit alignment with something that looks like the TML preamble
%                                 fprintf('----got TLM =%s ???\n',dec2bin(TLM_word,24));
                [ok,HOW_word]=paritycheck(obj.bit_sub_frame_buffer(1+30*(2-1):32+30*(2-1)));
%                                 fprintf('----got HOW =%s ???\n',dec2bin(HOW_word,24));
                D29D30=sum(obj.bit_sub_frame_buffer(32+30*(2-1)-1:32+30*(2-1)));
                if((ok)&&((D29D30==0)||(D29D30==2)))%for HOW the last 2 parity bits need to both be the same. this sign of them say if the signal is inverted or not. if 1,1 then we have an inverted signal
                    HOW_subframe_ID=bitand(bitshift(HOW_word,-2),7);
                    %yup matlab 2017a doesn't have hex literals
                    %TOW is for next subframe whitch is now as we have
                    %buffered 1 subrame and are on the start of the prn of the
                    %first bit of the next one
                    HOW_TOW=bitand(bitshift(HOW_word,-7),131071)*6;%0x1FFFF
                    obj.tow=HOW_TOW;
                    HOW_TOW_str=datestr(datetime(HOW_TOW,'ConvertFrom','epochtime','Epoch','1980-01-06'),'ddd HH:MM:SS');
                    %                     fprintf('----got TLM and HOW\n');
                    if(D29D30==2)
                        fprintf('signal is inverted\n');
                    end
                    fprintf('----Subframe ID=%d TOW: %s gps\n',HOW_subframe_ID,HOW_TOW_str);
                    fprintf('----%s\n',dec2bin(TLM_word,24));
                    fprintf('----%s\n',dec2bin(HOW_word,24));
                    obj.subframe(1)=TLM_word;
                    obj.subframe(2)=HOW_word;
                    obj.subframe_bad_word_count=0;
                    for word_n=3:10
                        [ok,word]=paritycheck(obj.bit_sub_frame_buffer(1+30*(word_n-1):32+30*(word_n-1)));
                        obj.subframe(word_n)=word;
                        if(ok)
                            fprintf('----%s\n',dec2bin(word,24));
                        else
                            fprintf('----%s bad\n',dec2bin(word,24));
                            obj.subframe_bad_word_count=obj.subframe_bad_word_count+1;
                        end
                    end
                    if(obj.subframe_bad_word_count==0)
                        fprintf('----YAY a good subframe\n');
                    else
                        fprintf("----TLM and HOW look ok but if this is a subframe then there are bad words\n");
                    end
                end
            end
        end
        
    end

end