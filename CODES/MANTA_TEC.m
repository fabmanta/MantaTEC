%% FUNCTION TO CALCULATE TEC FROM RINEX2.11
% Author :Fabio Manta 
% Last update: April-2020

%% Main Loop for station
% It returns a stroucture .mat file for each station with respect to all sat 
function MANTA_TEC(rinexDirectory,sp3fileName,ionexFileName,outputDirectory,...
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
global sv;
staz = []; % stations list
elapse_time = [];
eof = 0;

% h = waitbar(0,'do not rush!...');
h = waitbar(0);

 for ii = 1:length(a) %station loop
tic;
    staz = [staz; a(ii).name];

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
        % Keep only GPS satellites (G)
        name_sat = list_sat(strmatch('G',list_sat));
        numb_sat = size(name_sat,1);
        
        % create empty structure
        for index_obs = 1:numOfObsTypes
            subnum= char(typesOfObs(index_obs));
            Observables.(subnum)=nan(size(T,1),numb_sat);
        end

        % Fill the Observation structure
        for index_obs = 1:numOfObsTypes
            for sat = 1:numb_sat
                index_sat = strmatch(name_sat(sat),output(:,1));
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
        if (size(sp3fileName,1) > 1)
            sp3FileDyrectory = [rinexDirectory,sp3fileName(ii,:)];
        else
            sp3FileDyrectory = [rinexDirectory,sp3fileName];
        end
%         [SAT, sp3_list] = int_sp3(sp3FileDyrectory,tInterval,T(end));
        [SAT, sp3_list] = int_sp3(sp3FileDyrectory,T');

        % check for sv inconsistencies
        name_sat_srt =char(name_sat);
        sv =[];
        for sat =1:numb_sat
            sv(sat,:) =sscanf(name_sat_srt(sat,:),'G%d');
        end

        sp3_sv = [];
        exc_sv = [];
%         for isv = 1:length(sv)
%           if (find(sp3_list == sv(isv)))
            [sp3_sv,~] = intersect(sp3_list, sv);
%             sp3_sv = [sp3_sv sv(isv)];
%           else
            exc_sv = [setdiff(sp3_list, sv);setdiff(sv,sp3_list)];
            exc_sv = sort(exc_sv);
%             disp(['WARNING: PRN' num2str(sv(isv)) ' not in sp3 file']);
%             exc_sv = [exc_sv sv(isv)];
%           end
%         end

        % read ephemerides file: get transmitter group delay
        if (size(ionexFileName,1) > 1)
            ionsxFileDyrectory = [rinexDirectory,ionexFileName(ii,:)];
        else
            ionsxFileDyrectory = [rinexDirectory,ionexFileName];
        end
        if isempty(ionexFileName)
           tgd = [sv,ones(size(sv,1),1)];
           disp(['Warning: NO ionex file available']);
        else
            tgd = get_tgd_ionex(ionsxFileDyrectory);
        end

        % correct data for receiver interfrequency bias (added 3-16-2007, thomas dautermann)
%         ifb=ifbFromstruct(data,tgd);
        ifb =0;
        % get SIP information and create sip structure (only for svs present in sp3 file)
        [SIP_mat, emf_mat] = get_sip_fabio(SAT,sp3_sv,alt_ion,apcoords,Tmin,Tmax);

        % compute RAW TECU and time derivative
        if any(exc_sv)
            [~,id_remove]= intersect(sv,exc_sv);
%             id_remove = find(sv ==exc_sv);
            sv(id_remove) = [];
            L1(:,id_remove) = [];L2(:,id_remove) = [];
            C1(:,id_remove) = [];
            P1(:,id_remove) = [];P2(:,id_remove) = [];
            
            [~,id_tgd_remove] = intersect(tgd(:,1), exc_sv);
            tgd(id_tgd_remove,:) = [];
        end
            [IEC] = compute_tec(L1,L2,C1,P1,P2,sv,tgd,ifb,emf_mat);
        
        %calculate the distance between IPP and epicenter         
        [sat_distances,~]= distance_sat_station_calculator_3v(lon_epi,lat_epi,SIP_mat,altitutePPI,sv);
        
        % compute VTECU
%         iec_v = vertical_tecu_calculator_2v(IEC,SIP_mat,alt_ion,sv);    
        % compute Filt TECU
        IEC(find(IEC == 0)) = nan;
%         iec_v(find(iec_v == 0)) = nan;
        
        
        %correct for large steps in the data
        %You can use a group delay shift on the output to extract the relevant portion without that spike
        
%         [iec_filt,~] = filter_tecu(IEC,lowFreq,highFreq,tInterval);
%         [iec_v_filt,~] = filter_tecu(iec_v,lowFreq,highFreq,tInterval);
         iec_filt = filter_tecu_3v(IEC,lowFreq,highFreq,tInterval);
%          iec_v_filt = filter_tecu_3v(iec_v,lowFreq,highFreq,tInterval);
         iec_v_filt = vertical_tecu_calculator_2v(iec_filt,SIP_mat,alt_ion,sv);    

        
        %create IEC structure
        TEC_mat = [];
        for jj = 1:length(sv)
    %         iec_tmp = [epochs(:,4:6) iec(:,jj) iec_v(:,jj) iec_filt(:,jj) sat_distances(:,jj)];
            iec_tmp = [ IEC(:,jj) iec_v_filt(:,jj) iec_filt(:,jj) sat_distances(:,jj)];
            field = ['PRN' num2str(sv(jj))];
            TEC_mat = setfield(TEC_mat,field,iec_tmp);
        end

        % elevation cut-off
        [SIP_mat,TEC_mat] = elevation_cutoff(SIP_mat,TEC_mat,ele_cutoff,sv);
        
        %merrge the TEC and the SIP in the same structure
        tmp = merge_structure(t,SIP_mat,TEC_mat,sv);

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
save(fullfile(outputDirectory,fileName),'staz','SITE_COORD');


clear tmp;
end