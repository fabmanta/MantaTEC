%% FUNCTION TO CALCULATE TEC FROM RINEX2.11
% Author :Fabio Manta 
% Last update: March-2021
% 

%   https://gssc.esa.int/navipedia/index.php/Main_Page
%   Satellite system identifier
%    G or blank : GPS
%    R          : GLONASS
%    S          : Geostationary signal payload
%    T          : Transit
%    S          : SBAS payload
%    E          : Galileo
%    C          : BeiDou
%    J          : QZSS
%    I          : IRNSS/NavIC

%% Main Loop for station
% It returns a stroucture .mat file for each station with respect to all sat 
function rinex211_processor_2v(rinexDirectory,sp3fileName,ionexFileName,outputDirectory,...
    Year,defaultTimeInterval,lon_epi,lat_epi,altitutePPI,alt_ion,Tmin,Tmax,SITE_COORD,highFreq,lowFreq,ele_cutoff)

%     staz_NAME.PRNxx(01)  = Time
%     staz_NAME.PRNxx(02)  = LAT
%     staz_NAME.PRNxx(03) = LON
%     staz_NAME.PRNxx(04) = ELE
%     staz_NAME.PRNxx(05) = AZIM
%     staz_NAME.PRNxx(06) = EFM
%     staz_NAME.PRNxx(07)  = TEC
%     staz_NAME.PRNxx(08)  = TEC_V_FILT
%     staz_NAME.PRNxx(09)  = TEC_FILT
%     staz_NAME.PRNxx(10)  = DISTANCE

% =========================================================================
a = dir(sprintf('%s*.%so',rinexDirectory,Year(end-1:end)));

global staz;
staz = []; % stations list
elapse_time = [];
eof = 0;
% http://www2.unb.ca/gge/Resources/GLONASSConstellationStatus.txt
k = [1,-4,5,6,1,-4,5,6,-2,-7,0,-1,-2,-7,0,-1,4,-3,3,2,4,-3,3,2]; %GLONASS Constellation Status (Channels)

% h = waitbar(0,'do not rush!...');
h = waitbar(0);

 for ii = 1:length(a) %station loop
tic;
    staz = [staz; a(ii).name];
    tmp_k = k;
    
    if any(elapse_time)
        countdown = mean(elapse_time) * (length(a)-ii);
        countdown = countdown/3600;
        DateNumber = datenum(0,0,0,countdown,0,0);
        remainig_time = datestr(DateNumber,'HH:MM:SS');
        waitbar(ii/length(a),h,sprintf('Now procesing station %s %d/%d time remaining = %s',upper(staz(ii,1:4)),ii,length(a),remainig_time))
    else
        waitbar(ii/length(a),h,sprintf('Now procesing station %s %d/%d time remaining = . . .',upper(staz(ii,1:4)),ii,length(a)))
    end
    
    %   read rinex data file    
    %NEW version to read RINEX 2.11 just with matlab
    %read the header only
    [fid, rinexHeader, gnssType, markerName, anttInterval, apcoords,...
    numOfObsTypes, typesOfObs, tFirstObs, tLastObs, tInterval, ...
    timeSystem, numHeaderLines, clockOffsetsON, leapSec, eof] = ...
    rinexReadsObsFileHeader211(fullfile(rinexDirectory,staz(ii,:)));
    time_block = [];output = [];
    
    % check if tInterval exist
    if isnan(tInterval)
        tInterval = defaultTimeInterval;
    end
    
    %concatenate Rinex Obs Block to sngle matrix
    while 1
        line = fgetl(fid); % returns -1 if only reads EOF
        if line == -1
            break
        end
        %read the Block header only to retrive for each epoch: name and
        %number of the operative satellite.
        [fid,epochflag,clockOffset,date,numSV,SVlist] =rinexReadsObsBlockHead211(line,fid);
        [Obs, LLI, SS] =rinexReadsObsBlock211(fid, numSV, numOfObsTypes);


        time_block = [time_block; date(4) + date(5)/60 + round(date(6))/3600 *ones(numSV,1)];
        
        Obs_cel = num2cell(Obs);
        Obs_block = [ SVlist Obs_cel];
        output = [output;Obs_block];

    end
    fclose(fid);
    
    % Restructure the observation files
    % compute time vector (hours of current day)
    T = unique(time_block);
        
    % Check if the observation time is captured by the selected time window
    I = find(T>=Tmin & T<=Tmax);    
    if any(I)

        % find number of satellite
        list_sat = unique(output(:,1));
        type_sat= char(list_sat);
        type_sat= unique(type_sat(:,1));
        name_sat_srt = char(list_sat);

        % keep only GPS and GlONASS
        id_GPS = strmatch('G',type_sat);
        id_GLONASS = strmatch('R',type_sat); 
        if isempty(id_GLONASS)
            id_sat =  [strmatch(type_sat(id_GPS),name_sat_srt)];
        else
            id_sat =  [strmatch(type_sat(id_GPS),name_sat_srt);strmatch(type_sat(id_GLONASS),name_sat_srt)];
        end
        list_sat = list_sat(id_sat);
        name_sat_srt=name_sat_srt(id_sat,:);

        
        name_satG = list_sat(strmatch(type_sat(id_GPS),list_sat));
        name_satR = list_sat(strmatch(type_sat(id_GLONASS),list_sat));
        numb_sat = size(list_sat,1);
        numb_sat_G = size(name_satG,1);
        numb_sat_R = size(name_satR,1);

        % create empty structure
        for index_obs = 1:numOfObsTypes
            subnum= char(typesOfObs(index_obs));
            Observables.(subnum)=nan(size(T,1),numb_sat);
        end

        % Fill the Observation structure
        for index_obs = 1:numOfObsTypes
            for sat = 1:numb_sat
                index_sat = strmatch(list_sat(sat),output(:,1));
                [~,index_time] = ismember(T,time_block(index_sat)); 
                index_time = find(index_time);
                Observables.(typesOfObs{index_obs})(index_time,sat) = cell2mat(output(index_sat,index_obs+1));
            end
        end

        % fix length dataset filling time gaps
        Time_original = T;
        gappen_original = find(diff(T) > (tInterval/3600)*1.1);
        if any(gappen_original)
            gap_sec = (T(gappen_original+1)-T(gappen_original))*3600;
            n_missing_points = round((gap_sec/tInterval)-1);
            gappen = gappen_original;
            new_gaps = [];
%             while ~isempty(gappen)
            for gg = 1:size(gappen,1)
                new_gaps = [new_gaps;gappen(1)];
                A = T(1:gappen(1));
                B = T(gappen(1)+1:end); 
                T = [A;T(gappen(1))+(tInterval*(1:n_missing_points(gg))/3600)';B];
                gappen = find(diff(T) > (tInterval/3600)*1.1);
            end   
            tmp = Observables;
            % create empty structure
            for index_obs = 1:numOfObsTypes
                subnum= char(typesOfObs(index_obs));
                Observables.(subnum)=nan(size(T,1),numb_sat);
            end

            [~,index_time] = ismember(T,Time_original);
            index_time = find(index_time);
            for index_obs = 1:numOfObsTypes
                for sat = 1:numb_sat    
                Observables.(typesOfObs{index_obs})(index_time,sat) ...
                    = tmp.(typesOfObs{index_obs})(:,sat);       
                end
            end
        end

        % Organize observables into separate arrays
        L1=Observables.L1; L2=Observables.L2;
        P2=Observables.P2; 
        if (isfield(Observables,'P1')) P1=Observables.P1; end
        if (isfield(Observables,'C1')) C1=Observables.C1; end

        % replace NaN by zeros
        L1(isnan(L1)) = 0; L2(isnan(L2)) = 0;
        P2(isnan(P2)) = 0; 

        if (isfield(Observables,'P1'))
           P1(isnan(P1)) = 0;
        else
           P1 = [];
        end

        if (isfield(Observables,'C1'))
            C1(isnan(C1)) = 0;
        else
            C1 = zeros(size(L1));
        end

        % if P1 not defined in rinex or P1 defined but empty, then no P1 data => use C1 instead
        if (isempty(find(P1)) || length(find(P1)) < length(find(C1)))
          P1 = C1;
        %       disp(['WARNING: site ' nm_site ' does not have P1 data --> using C1']);
        end

        % select time interval
        I = find(T>=Tmin & T<=Tmax);
        t = T(I); 
        L1 = L1(I,:); L2 = L2(I,:); 
        C1 = C1(I,:);
        P1 = P1(I,:); P2 = P2(I,:);     

        %filter pr and phase oscillations
        clean_data = 'y';
        if clean_data=='y'
             disp('Cleaning Observables')
         [L1,L2,C1,P1,P2]=clean_obs(L1,L2,C1,P1,P2);
        end %if clean_data


        % get coordinates of GPS site
        W=xyz2wgs([0 apcoords']);
        lat_site = W(3); lon_site = W(2); alt_site = 0;
        SITE_COORD=[SITE_COORD; [W(3) W(2)]];

        % read sp3file and interpolate every splint seconds 
        if gnssType == 'M'
        for kk = 1:size(type_sat,1)
            switch type_sat(kk)
                case 'R' % GLONASS
                sp3FileDyrectory = [rinexDirectory,sp3fileName(1,:)];
                [SAT_R, sp3_listR] = int_sp3(sp3FileDyrectory,T');
                case 'G' % GPS
                sp3FileDyrectory = [rinexDirectory,sp3fileName(2,:)];
                [SAT_G, sp3_listG] = int_sp3(sp3FileDyrectory,T');
            end
        end
        else
            sp3FileDyrectory = [rinexDirectory,sp3fileName(2,:)];
            [SAT_G, sp3_listG] = int_sp3(sp3FileDyrectory,T');
        end
        
        % check for sv inconsistencies
        name_sat_srtG =char(name_satG);
        name_sat_srtR =char(name_satR);
        sv_G =[];
        sv_R =[];


        for sat_G =1:numb_sat_G
            sv_G(sat_G,:) =sscanf(name_sat_srtG(sat_G,:),strcat(type_sat(id_GPS),'%d'));
        end
        
        if ~isempty(id_GLONASS)
            for sat_R =1:numb_sat_R
                sv_R(sat_R,:) =sscanf(name_sat_srtR(sat_R,:),strcat(type_sat(id_GLONASS),'%d'));
            end
        end
        tot_sv= [sv_G;sv_R];

        % GPS
        sp3_G = [];
        exc_G = [];
        [sp3_G,~] = intersect(sp3_listG, sv_G);
        exc_G = [setdiff(sp3_listG, sv_G);setdiff(sv_G,sp3_listG)];
        exc_G = sort(exc_G);
        
        if ~isempty(id_GLONASS)
            % GLONASS
            sp3_R = [];
            exc_R = [];
            [sp3_R,~] = intersect(sp3_listR, sv_R);
            if ~isempty(sv_R)
                exc_R = [setdiff(sp3_listR, sv_R);setdiff(sv_R,sp3_listR)];
                exc_R = sort(exc_R);
            end
        end
        
        % read ephemerides file: get transmitter group delay
        if (size(ionexFileName,1) > 1)
            ionsxFileDyrectory = [rinexDirectory,ionexFileName(ii,:)];
        else
            ionsxFileDyrectory = [rinexDirectory,ionexFileName];
        end
        if isempty(ionexFileName)
           tgd_G = [sv_G,ones(size(sv_G,1),1)];
           disp(['Warning: NO ionex file available']);
        else
            tgd_G = get_tgd_ionex(ionsxFileDyrectory);
            tgd_R = [sv_R,ones(size(sv_R,1),1)];
        end

        % correct data for receiver interfrequency bias (added 3-16-2007, thomas dautermann)
%         ifb=ifbFromstruct(data,tgd);
        ifb =0;
        % get SIP information and create sip structure (only for svs present in sp3 file)
        [SIP_mat_G, emf_mat_G] = get_sip_fabio(SAT_G,sp3_G,alt_ion,apcoords,Tmin,Tmax);
        if ~isempty(id_GLONASS)
           [SIP_mat_R, emf_mat_R] = get_sip_fabio(SAT_R,sp3_R,alt_ion,apcoords,Tmin,Tmax);
        else
            
        end
        
        % compute RAW TECU and time derivative
        tot_exc =[exc_G;exc_R];
        if any(tot_exc)
            [~,id_removeG]= intersect(sv_G,exc_G);
            [~,id_removeR]= intersect(sv_R,exc_R);
            id_remove = [id_removeG;id_removeR+length(sp3_G)];
%             id_remove = find(sv ==exc_sv);
            list_sat(id_remove,:) = [];
            tot_sv(id_remove) = [];
            sv_G(id_removeG) = [];
            sv_R(id_removeR) = [];
            L1(:,id_remove) = [];L2(:,id_remove) = [];
            C1(:,id_remove) = [];
            P1(:,id_remove) = [];P2(:,id_remove) = [];

            [~,id_tgd_remove] = intersect(tgd_G(:,1), exc_G);
            tgd_G(id_tgd_remove,:) = [];
            [~,id_tgd_remove] = intersect(tgd_R(:,1), exc_R);
            tgd_R(id_tgd_remove,:) = [];
        end
        
        % remove exciding GLONASS satellite from k channel list;
        exciding = setdiff(1:24,sv_R);
        tmp_k(exciding) = [];

        [IEC_G] = compute_tec_2v(L1(:,1:length(sv_G)),L2(:,1:length(sv_G)),C1(:,1:length(sv_G)),P1(:,1:length(sv_G)),P2(:,1:length(sv_G)),sv_G,tgd_G,ifb,emf_mat_G,type_sat(id_GPS),0);
        IEC_G(find(IEC_G == 0)) = nan;        
        %calculate the distance between IPP and epicenter         
        [sat_distances_G,~]= distance_sat_station_calculator_3v(lon_epi,lat_epi,SIP_mat_G,altitutePPI,sv_G);
        %filt stec and vtec
        iec_filt_G = filter_tecu_3v(IEC_G,lowFreq,highFreq,tInterval);
        iec_v_filt_G = vertical_tecu_calculator_2v(iec_filt_G,SIP_mat_G,alt_ion,sv_G);    
         
        if (gnssType == 'M') && ~isempty(id_GLONASS)

            [IEC_R] = compute_tec_2v(L1(:,length(sv_G)+1:end),L2(:,length(sv_G)+1:end),C1(:,length(sv_G)+1:end),P1(:,length(sv_G)+1:end),P2(:,length(sv_G)+1:end),sv_R,tgd_R,ifb,emf_mat_R,type_sat(id_GLONASS),tmp_k);
            IEC_R(find(IEC_R == 0)) = nan;
            [sat_distances_R,~]= distance_sat_station_calculator_3v(lon_epi,lat_epi,SIP_mat_R,altitutePPI,sv_R);
            iec_filt_R = filter_tecu_3v(IEC_R,lowFreq,highFreq,tInterval);
            iec_v_filt_R = vertical_tecu_calculator_2v(iec_filt_R,SIP_mat_R,alt_ion,sv_R);    
        else
          IEC_R =[];  
          sat_distances_R=[];
          iec_filt_R=[];
          iec_v_filt_R=[];
        end  

        %create TEC structure
        TEC_mat = [];
        SIP_mat = [];
        for jj = 1:length(sv_G)
            iec_tmp = [ IEC_G(:,jj) iec_v_filt_G(:,jj) iec_filt_G(:,jj) sat_distances_G(:,jj)];
            field = sprintf('G%.2d', sv_G(jj));
            TEC_mat = setfield(TEC_mat,field,iec_tmp);
            SIP_mat = setfield(SIP_mat,field,eval(sprintf('SIP_mat_G.PRN%d',sv_G(jj))) );
        end
        for jj = 1:length(sv_R)
            iec_tmp = [ IEC_R(:,jj) iec_v_filt_R(:,jj) iec_filt_R(:,jj) sat_distances_R(:,jj)];
            field = sprintf('R%.2d', sv_R(jj));
            TEC_mat = setfield(TEC_mat,field,iec_tmp);
            SIP_mat = setfield(SIP_mat,field,eval(sprintf('SIP_mat_R.PRN%d',sv_R(jj))) );
        end

        % elevation cut-off
        [SIP_mat,TEC_mat] = elevation_cutoff_2v(SIP_mat,TEC_mat,ele_cutoff,name_sat_srt,tot_sv);


        % merge the TEC and the SIP in the same structure
        tmp = merge_structure_2v(t,SIP_mat,TEC_mat,char(list_sat));

        v = genvarname(eval('sprintf(staz(end,1:end-4))'));
        eval([v ' = tmp']);

        %Save all
            fileName = sprintf('%s',staz(ii,1:end-4));

            if exist(outputDirectory,'dir')
                save(fullfile(outputDirectory,fileName),sprintf('%s',staz(ii,1:end-4)));
            else
                mkdir(outputDirectory)
                save(fullfile(outputDirectory,fileName),sprintf('%s',staz(ii,1:end-4)));
            end

         
    else
        disp(['WARNING: station ' upper(staz(ii,1:4)) ' is not recording at the time of the event']);

    end
        
clear Observables;

elapse_time = [elapse_time, toc];
 end
 
close(h)


fileName = 'staz_coordinates';
name_sat_srt = char(list_sat);
save(fullfile(outputDirectory,fileName),'staz','SITE_COORD','name_sat_srt');


clear tmp;
end