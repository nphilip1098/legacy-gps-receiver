classdef satposcalulator_class < handle
    %calculate sat positions
    %

    properties
        ECEF;
        sv_clock_bias;
    end
    
    properties (GetAccess = public , SetAccess = private)
        
        ephemerides;
        
        sv_ephemeris;
        
        svid_;
        m0_;
        dn_;
        e_ ;
        a_ ;
        omg0_;
        i0_ ;
        w_ ;
        odot_;
        idot_;
        cuc_;
        cus_;
        crc_;
        crs_;
        cic_;
        cis_;
        toe_;
        iode_;
        GPS_week_;
        toc_;
        af0_;
        af1_;
        af2_;
        TGD_;
        
    end
    
    methods
        
        function read_nav_file(obj,filename)
            % reads a file or a set of rinex 2.x gps nav files
            if(iscell(filename))
                obj.ephemerides = read_rinex_nav(char(filename(1)));
                for k=2:numel(filename)
                    obj.ephemerides = [obj.ephemerides;read_rinex_nav(char(filename(k)))];
                end
            else
                obj.ephemerides = read_rinex_nav(filename);
            end
        end
        
        function Reset(obj)
            obj.ECEF=[0,0,0];
            obj.sv_clock_bias=0;
        end
    
        function obj = satposcalulator_class()
            obj.Reset();
        end
        
        %towreq is sv time received for L1 receivers if correct_for_sv_time==true else it's gps time
        function [ ECEF ] = getpos(obj,sv,weekreq,towreq,correct_for_sv_time )
            %calculates position of sv given tsv,ephemeris,towreq return ECDF position
            % Input values:
            % sv: sv number
            % weekreq: week of the request in GPS time
            % towreq: time of week request in GPS time
            % correct_for_sv_time: if true then will correct towrequest for relitivity error and sv clock bias
    		% NB all ephemerides need to be in the same week as the request
            
            %get the closest ephemeride wrt weekreq and towreq.
            %will have week rollover, I think in 2080 due to cal2gpstime
            sv_ephemerides = obj.ephemerides(obj.ephemerides(:,1)==sv,:);
            [~,sv_ephemeride_index]=min(abs((towreq-sv_ephemerides(:,17))+(604800*(weekreq-sv_ephemerides(:,19)))));
            obj.sv_ephemeris=sv_ephemerides(sv_ephemeride_index,:);
            toedelta=(towreq-obj.sv_ephemeris(17));
            assert(abs(toedelta)<(4*60*60),'toe wrt now has to be less than 4 hours');
            
            %just the best ephemeris for this sv
            ephemeris=obj.sv_ephemeris;
            
            GM = 3.986005e14;             % earth's universal gravitational [m^3/s^2]
            c = 2.99792458e8;             % speed of light (m/s)
            omegae_dot = 7.2921151467e-5; % earth's rotation rate (rad/sec)
            
            % initialize constants and variables
            svid = ephemeris(:,1);
            m0   = ephemeris(:,2);
            dn   = ephemeris(:,3);
            e    = ephemeris(:,4);
            a    = (ephemeris(:,5)).^2;
            omg0 = ephemeris(:,6);
            i0   = ephemeris(:,7);
            w    = ephemeris(:,8);
            odot = ephemeris(:,9);
            idot = ephemeris(:,10);
            cuc  = ephemeris(:,11);
            cus  = ephemeris(:,12);
            crc  = ephemeris(:,13);
            crs  = ephemeris(:,14);
            cic  = ephemeris(:,15);
            cis  = ephemeris(:,16);
            toe  = ephemeris(:,17);
            iode = ephemeris(:,18);
            GPS_week = ephemeris(:,19);
            toc=ephemeris(:,20);
            af0= ephemeris(:,21);
            af1= ephemeris(:,22);
            af2= ephemeris(:,23);
            TGD=ephemeris(:,24);
            
            %done this because matlab doesn't do code highlighting on
            %properites
            obj.svid_=svid;
            obj.m0_=m0;
            obj.dn_=dn;
            obj.e_=e;
            obj.a_=a;
            obj.omg0_=omg0;
            obj.i0_=i0;
            obj.w_=w;
            obj.odot_=odot;
            obj.idot_=idot;
            obj.cuc_=cuc;
            obj.cus_=cus;
            obj.crc_=crc;
            obj.crs_=crs;
            obj.cic_=cic;
            obj.cis_=cis;
            obj.toe_=toe;
            obj.iode_=iode;
            obj.GPS_week_=GPS_week;
            obj.toc_=toc;
            obj.af0_=af0;
            obj.af1_=af1;
            obj.af2_=af2;
            obj.TGD_=TGD;
            
            nn = size(ephemeris,1);
            for ii = 1:nn
                
                %deal with week rollover.
                if((towreq-toe(ii))>(604800/2))
                    towreq=towreq-604800;
                elseif((toe(ii)-towreq)>(604800/2))
                    towreq=towreq+604800;
                end
                
                % Procedure for coordinate calculation
                n0 = sqrt(GM/a(ii)^3); % (rad/s)
                tk = towreq-toe(ii);      % Time from eph ref epoch (s)
                n = n0+dn(ii);         % Corrected mean motion (rad/s)
                M = m0(ii)+n*tk;       % Mean anomaly (rad/s)
                
                obj.sv_clock_bias=0;
                if(correct_for_sv_time)
                    % Perform Newton-Raphson solution for Eccentric anomaly estimate (rad)
                    NRnext = 0;
                    NR = 1;
                    m = 1;
                    while (abs(NRnext-NR)>1e-15)
                        NR=NRnext;
                        f=NR-e(ii)*sin(NR)-M;
                        f1=1-e(ii)*cos(NR);
                        f2=e(ii)*sin(NR);
                        NRnext=NR-(f/(f1-(f2*f/2*f1)));
                        m=m+1;
                    end
                    E=NRnext; % Eccentric anomaly estimate for computing delta_tr (rad)
                    
                    % Time correction
                    F = -2*sqrt(GM)/c^2; % (s/m^1/2)
                    delta_tr = F*e(ii)*sqrt(a(ii))*sin(E);% relativity eccentricity error
                    delta_tsv = af0(ii)+af1(ii)*(towreq-toe(ii))+delta_tr;%sv clock offset
                    delta_tsv = delta_tsv -TGD;%jonti for L1C/A
                    t = towreq-delta_tsv;
                    obj.sv_clock_bias=delta_tsv;
                    
                    %deal with week rollover.
                    if((t-toe(ii))>(604800/2))
                        t=t-604800;
                    elseif((toe(ii)-t)>(604800/2))
                        t=t+604800;
                    end
                    
                    tk=t-toe(ii);  		 % Time from eph ref epoch (s)
                    M=m0(ii)+n*tk;	     % Mean anomaly (rad/s)
                end
                
                % Perform Newton-Raphson solution for Eccentric anomaly (rad)
                NRnext=0;
                NR=1;
                m=1;
                while (abs(NRnext-NR)>1e-15)
                    NR=NRnext;
                    f=NR-e(ii)*sin(NR)-M;
                    f1=1-e(ii)*cos(NR);
                    f2=e(ii)*sin(NR);
                    NRnext=NR-(f/(f1-(f2*f/2*f1)));
                    m=m+1;
                end
                E=NRnext; % Eccentric anomaly (rad)
                v = atan2(sqrt(1-e(ii)^2)*sin(E), cos(E)-e(ii));
                phi = v+w(ii);
                u = phi                    + cuc(ii)*cos(2*phi)+cus(ii)*sin(2*phi);
                r = a(ii)*(1-e(ii)*cos(E)) + crc(ii)*cos(2*phi)+crs(ii)*sin(2*phi);
                i = i0(ii)+idot(ii)*tk     + cic(ii)*cos(2*phi)+cis(ii)*sin(2*phi);
                Omega = omg0(ii)+(odot(ii)-omegae_dot)*tk-omegae_dot*toe(ii);
                x1 = cos(u)*r;
                y1 = sin(u)*r;
                
                % ECEF coordinates
                assert(svid(ii)==sv,"error sv doesn't match");
                ECEF(1,ii) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
                ECEF(2,ii) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
                ECEF(3,ii) = y1*sin(i);
                
            end
            
            obj.ECEF=ECEF;
            
        end
        
    end

end