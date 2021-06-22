%% MAIN SCRIPT TO CALCULATE TEC FROM RINEX FILES
% Author :Fabio Manta 
% Last update: June-2021

% Data centers GNSS 
% https://www.unavco.org/data/gps-gnss/data-access-methods/dai2/app/dai2.html#
% ftp://cddis.gsfc.nasa.gov/gnss/data/daily/
% http://stdb2.isee.nagoya-u.ac.jp/GPS/GPS-TEC/gnss_provider_list.html
% ftp://gssc.esa.int/gnss/data/daily/
% ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/
% ftp://nfs.kasi.re.kr/gps/data/
% ftp://igs.gnsswhu.cn/pub/gnss/obs/
% ftp://ftp.earthobservatory.sg/SugarData
% ftp://gpsfree.gm.ingv.it/OUTGOING/RINEX30/RING/
% http://garner.ucsd.edu/pub/rinex/
% ftp://www.ordnancesurvey.co.uk/gps/rinex
% ftp://rinex.ordnancesurvey.co.uk/gps/rinex

clear all;close all;clc;
warning off

% =========================================================================
% Set up Directories and Event info

global staz; % list of GNSS stations
defaultTimeInterval = 30; %aquisition period in seconds
% Frequency range for the Butterworth filter
lowFreq = 0.001;% Hz
highFreq = 0.010;%Hz; 
ele_cutoff = 15; % Elevation cut-off in degrees from the horizon

evntName = 'Mw8.3'; %'sinabung' or 'Mw7.0'
dateString = '20150916';
lat_epi = -31.57; %Latitude epicenter;
lon_epi = -71.654; %Longitude epicenter;
altitutePPI = 300; %ionosphere altitude in km
alt_ion = altitutePPI *1000; %ionosphere altitude in meters
Tmin= 22; Tmax= 24;% Time range to be analyze in hours (UTC time)
SITE_COORD = []; %coordinates of the GNSS stations 
EventTime = datenum(2015,09,16,22,54,32); % UTC time 
[Year,Month,Day,Hours,Minutes,Seconds]= datevec(EventTime);
Year = num2str(Year); 
EventTimeHours = Hours+Minutes/60+Seconds/3600;
% The following variables are necessary to set the limits of the xaxis in
% the plots
startTime = round(EventTimeHours-1);% in UTC time
endTime = 24;%

root_direct = '/Users/fabio/Dropbox/MANTA-TEC_toolkit';
addpath(fullfile(root_direct,'/codes'))
addpath(root_direct)

rinexFolder = sprintf('%s_%s_%dkm',dateString,evntName,altitutePPI);
outputFolder = sprintf('output_%s_%s_%dkm',dateString,evntName,altitutePPI);
plotFolder = sprintf('plots_%s_%s_%dkm',dateString,evntName,altitutePPI);

rinexDirectory = sprintf('%s/input/%s/',root_direct,rinexFolder);
outputDirectory = sprintf('%s/output/%s/',root_direct,outputFolder);
plotsDirectory = sprintf('%s/plots/%s/',root_direct,plotFolder);

% http://www.gnsscalendar.com/index.html?year=2018
% SP3 files are in GPS week
sp3 = dir(sprintf('%s*.sp3',rinexDirectory));
sp3 = struct2cell(sp3);
sp3fileName = cell2mat(sp3(1,:)'); %ftp://igs.ensg.ign.fr/pub/igs/products/

ionex = dir(sprintf('%s*.%si',rinexDirectory,Year(end-1:end)));
ionex = struct2cell(ionex);
ionexFileName = cell2mat(ionex(1,:)'); %ftp://ftp.aiub.unibe.ch/CODE/
% sp3FileDyrectory = [rinexDirectory,sp3fileName];
% ionsxFileDyrectory = [rinexDirectory,ionexFileName];

spectro = 'spectro/';
spectroDirectory = [plotsDirectory,spectro];
mkdir(outputDirectory)
mkdir(spectroDirectory)

%% Calculate all parameters 
% It returns a structures .mat file for each station with respect to all sat 

%     nameDOY.PRNxx(:,01)  = Time
%     nameDOY.PRNxx(:,02)  = LAT_IPP
%     nameDOY.PRNxx(:,03) = LON_IPP
%     nameDOY.PRNxx(:,04) = ELE
%     nameDOY.PRNxx(:,05) = AZIM
%     nameDOY.PRNxx(:,06) = EFM
%     nameDOY.PRNxx(:,07) = STEC
%     nameDOY.PRNxx(:,08) = VTEC_FILT
%     nameDOY.PRNxx(:,09) = STEC_FILT
%     nameDOY.PRNxx(:,10) = DISTANCE 

MANTA_TEC(rinexDirectory,sp3fileName,ionexFileName,outputDirectory,...
                   Year,defaultTimeInterval,lon_epi,lat_epi,altitutePPI,alt_ion,Tmin,Tmax,SITE_COORD,highFreq,lowFreq,ele_cutoff);

% Load data of the Day of the event
cd(outputDirectory)
load_file_tec
load('staz_coordinates.mat')
cd(plotsDirectory)
% =========================================================================
%% DATA VISUALIZATION
%% HODOCHRONES

TECU_delta = 0.5;
save_fig = 'y'; % possyble values 'y' or 'n'
SatNum = (1:32); % all sats;
% SatNum = [31];

for sat_loop = 1:length(SatNum) %loop satellites
    
    figure
    set(gcf,'color','w');
 
    for staz_loop = 1: size(staz,1) %loop stations
            
        if (isfield(eval(staz((staz_loop),1:end-4)),sprintf('PRN%d',SatNum(sat_loop))))

            tmpTIME = eval(strcat(staz((staz_loop),1:end-4),sprintf('.PRN%d(:,1)',SatNum(sat_loop)))); % related with Lat Lon
            tmpDistance = eval(strcat(staz((staz_loop),1:end-4),sprintf('.PRN%d(:,10)',SatNum(sat_loop))));
            tmpTecu = eval(strcat(staz((staz_loop),1:end-4),sprintf('.PRN%d(:,8)',SatNum(sat_loop))));                       
%             tmpTecu = signal_attenuation(tmpTecu);
%             [SIP_mat,TEC_mat] = elevation_cutoff(SIP_mat,TEC_mat,ele_cutoff,sv);

            SatTrajectory(staz_loop) = cplot(tmpTIME,tmpDistance,tmpTecu,'o','markerfacecolor','flat','MarkerSize',5); %'linewidth',5);%
            hold on 
       
        end      
%             pause()
    end
        plot([EventTimeHours, EventTimeHours],[0, 1500],'--','color',[.7 .7 .7]);
        plot([EventTimeHours, EventTimeHours+0.5556],[0, 2000],'--','color',[1 0 0]);   %1 km/s

        set(gca, 'CLim', [-TECU_delta,TECU_delta]);
        set(gca,'FontSize',20)
        colorbar;
        colormap(jet)
        
        set(gca,'XTick',(startTime:1:endTime))
        axis([startTime endTime 0 1000])
        xlabel('Time','FontSize', 20);
        ylabel('Epicentral Range (Km)','FontSize', 20);
        title(sprintf('VTEC-G%d',SatNum(sat_loop)))
        
end 

% SAVE HODOCHRONES

        if (save_fig == 'y')
                for n = 1:length(SatNum)
                    namesatellite = SatNum(n);
                    saveas(figure(n),sprintf('%sAll_staz_G%d',plotsDirectory,namesatellite),'fig')
                    print(figure(n),'-depsc','-r0',sprintf('%sAll_staz_G%d',plotsDirectory,namesatellite));
                    print(gcf,'-dpng','-r0',sprintf('%sAll_staz_G%d',plotsDirectory,namesatellite));
                end
        close all
        end
        
%% FIND CLOSER STATIONS TO SOURCE

% conversion to radian epicentral coordinates
lon_epi_rad=lon_epi*pi/180;
lat_epi_rad=lat_epi*pi/180;
radius=6371; %in Km

% loop for all the stations
for jj = 1:length(staz)
    % get stations coordinnates and covert them in radian
    lat_staz = SITE_COORD(jj,1);
    lon_staz = SITE_COORD(jj,2);
    lon_staz_rad=lon_staz*pi/180;
    lat_staz_rad=lat_staz*pi/180;

    xx=acos(sin(lat_epi_rad).*sin(lat_staz_rad)+cos(lat_epi_rad).*cos(lat_staz_rad).*cos((lon_staz_rad-lon_epi_rad)));
    dist_staz_epi(:,jj) = radius.*xx;% calculate GNSS distance to the epicenter
end

id_st = find(dist_staz_epi <= 500); % id of the stations at less than Xkm from the epicenter
id_lat =  find(SITE_COORD(:,1) <= lat_epi + 1); % id of the stations sother then the epicenter 

[val,pos]=intersect(id_st,id_lat); %keep only only the station id that satisfy both above criteria
        
%% TRAJECTORY MAP
%Ionospheric Pierce Point trajectory

id_staz= [1:size(staz,1)]; % all GNSS id selected
% id_staz= [2]; % select only specific GNSS id 
% id_staz = val'; % select only the GNSS that satisfy the criteria defined in the previous section
% SatNum = [5 12 13 15 17 19 20]; % here you can re-define the set of satellites to plot in the map

% Set the map extent
delta_deg = 20; 
lat_ex=[lat_epi-10 lat_epi+delta_deg]; %
lon_ex=[lon_epi-delta_deg lon_epi+delta_deg];

m_proj('Mercator','longtitudes',lon_ex,'latitudes',lat_ex);
m_grid('linest','none','box','fancy','tickdir','in','fontsize',16);
%for higher resolution use gshhs_i, for full resolution use gshhs_f
m_gshhs_h('patch',[1.0000    1.0000    0.8750]);

hold on
[xEvent, yEvent]= m_ll2xy(lon_epi,lat_epi);
plot(xEvent,yEvent,'pr','MarkerSize',20,'MarkerFaceColor', 'r');
[xstaz, ystaz]= m_ll2xy(SITE_COORD(:,2),SITE_COORD(:,1));
plot(xstaz(id_staz),ystaz(id_staz),'ok','MarkerSize',4,'MarkerFaceColor', 'k'); 
text(xstaz(id_staz)+0.001,ystaz(id_staz)+0.001,upper(staz(id_staz,1:4)), 'FontSize', 10);

headWidth = 2;
headLength = 2;

marker = [ 'o' 's' 'd' 'v'];
colore = ['m' 'g' 'b' 'c'];
set(gca,'FontSize',16)
set(gcf,'color','w');
box on  
  
 for jj= 1: size(id_staz,2) 
         
        for ii = 1:length(SatNum) 
            
          hold on    

            if (isfield(eval(staz(id_staz(jj),1:end-4)),sprintf('PRN%d',SatNum(ii))))

            tmpTIME = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,1)',SatNum(ii)))); % related with Lat Lon
            tmpLAT = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,2)',SatNum(ii))));  
            tmpLON = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,3)',SatNum(ii))));
            [tmpLON, tmpLAT] = m_ll2xy(tmpLON,tmpLAT);

            startingPosition = find(tmpTIME <= EventTimeHours-0.5);
            finalPosition = find(tmpTIME <= EventTimeHours+1);
            eventPosition = find(tmpTIME <= EventTimeHours);

          
                    if ~isempty (startingPosition) 
                        if startingPosition(end)==length(tmpLON)                        
                             text(tmpLON(end),tmpLAT(end),sprintf('G%d',SatNum(ii)));
                        else
                            startingPosition = startingPosition(end);
                            finalPosition = finalPosition(end);
                            eventPosition = eventPosition(end);
                            
                            plot(tmpLON(startingPosition:finalPosition),tmpLAT(startingPosition:finalPosition),'Color',[.7 .7 .7])
                            if (size(SatNum,2)<=4)
                                cc = [colore(ii),marker(ii)];
                                plot(tmpLON(eventPosition),tmpLAT(eventPosition),cc,'MarkerFaceColor',colore(ii),'MarkerEdgeColor','k','MarkerSize',9)
                                text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('G%d',SatNum(ii)));
                            else
                              plot(tmpLON(eventPosition),tmpLAT(eventPosition),'ko','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','k','MarkerSize',9)
                              text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('G%d',SatNum(ii)));
                            end
                            ah = annotation('arrow',...
                                'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth,'Color',[.7 .7 .7]);
                             set(ah,'parent',gca);
                             if ~isnan(tmpLON(finalPosition))
                             set(ah,'position',[tmpLON(finalPosition) tmpLAT(finalPosition) tmpLON(finalPosition)-tmpLON(finalPosition-1) tmpLAT(finalPosition)-tmpLAT(finalPosition-1)]);
                             end
                        end
                    else
                        startingPosition = 1;
                        finalPosition = finalPosition(end);
                        eventPosition = eventPosition(end);
                        plot(tmpLON(startingPosition:finalPosition),tmpLAT(startingPosition:finalPosition),'Color',[.7 .7 .7])
                        if (size(SatNum,2)<=4 )
                            cc = [colore(ii),marker(ii)];
                            plot(tmpLON(eventPosition),tmpLAT(eventPosition),cc,'MarkerFaceColor',colore(ii),'MarkerEdgeColor','k','MarkerSize',9)
                            text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('G%d',SatNum(ii)));
                        else
                            cc = [colore(1),marker(1)];    
                            plot(tmpLON(eventPosition),tmpLAT(eventPosition),cc,'MarkerFaceColor',colore(1),'MarkerEdgeColor','k','MarkerSize',9)
                            text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('G%d',SatNum(ii)));                    
                        end
%                             text(tmpLON(finalPosition),tmpLAT(finalPosition),sprintf('G%d',SatNum(ii)+1),'FontSize',6);
                        ah = annotation('arrow',...
                            'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth,'Color',[.7 .7 .7]);
                         set(ah,'parent',gca);
                         if ~isnan(tmpLON(finalPosition-1))
                         set(ah,'position',[tmpLON(finalPosition) tmpLAT(finalPosition) tmpLON(finalPosition)-tmpLON(finalPosition-1) tmpLAT(finalPosition)-tmpLAT(finalPosition-1)]);
                         end                        
                    end
            end
       end

 end          

daspect([1 1 1])
        
saveas(gcf,'map','fig')
print(gcf,'-depsc','-r0','map');
print(gcf,'-dpng','-r0','map');
                       
%% FILT TEC timeseries
% Only stations paired with satellite N in order of distance form the epicenter at the time of the event

id_selection= id_staz; % The user can re-define the ID of the set of sations to be plot
staz_new = staz(id_selection,:);
yLowerLim = -0.4; % set Y axis lower limit

for ii = 1:size(staz_new,1) % loop for the selected stations
     dist_staz =[];

    for jj = 1:length(SatNum) % loop for the selected satellites
         
        % only if the satellite exist the retrive the value of IPP-epicenter distance
        if (isfield(eval(staz_new((ii),1:end-4)),sprintf('PRN%d',SatNum(jj))))

            tmpTIME = eval(strcat(staz_new((ii),1:end-4),sprintf('.PRN%d(:,1)',SatNum(jj)))); 
            tmpDistance = eval(strcat(staz_new((ii),1:end-4),sprintf('.PRN%d(:,10)',SatNum(jj)))); 
            % check if at the time of the event the satellite was close to the source
            eventPosition = find(tmpTIME <= EventTimeHours);
            if isempty(eventPosition)
                eventPosition = 1;
            else 
                eventPosition = eventPosition(end);
            end
            
            % if the satellite was close to the source and recording, then
            % save the value of IPP-epicenter distance
            dist_mat(ii,jj) = round(tmpDistance(eventPosition)); 
        else
            dist_mat(ii,jj) = NaN;
        end

    end
end

numberSatellites = size(SatNum,2);
distVect = reshape(dist_mat.',1,[]);
[value, id_sort] = sort(distVect); % sort the GNSS-IPP for the values of distance from the epicenter
ind = ceil(id_sort/numberSatellites);

figure
a=0;
ticks = [];        
staz_name = [];

for jj = 1:length(ind)
    
    id_sat = id_sort(jj) -(ind(jj)-1)*numberSatellites;
    
    % PLOT ONLY THE GNSS-SAT PAIRS WHICH WERE RECORDING AT THE TIME OF THE EVENT
    if (isfield(eval(staz_new(ind(jj),1:end-4)),sprintf('PRN%d',SatNum(id_sat)))) && ~isnan(value(jj))

        tmpTIME = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.PRN%d(:,1)',SatNum(id_sat)))); 
        tmpTEC = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.PRN%d(:,8)',SatNum(id_sat)))); 
        tmpTEC = signal_attenuation(tmpTEC);
        hold on;

        ticks = [ticks, a];
        pairs = [staz_new(ind(jj),1:4),sprintf('-PRN%2d',SatNum(id_sat))];
        staz_name = [staz_name; pairs];
        plot(tmpTIME,tmpTEC+a,'r')
        a = a + abs(yLowerLim);
    end

end

yyaxis left
ylimt = get(gca,'Ylim');
plot([EventTimeHours, EventTimeHours],[ylimt(1), ylimt(2)],'--','color',[0.1 0.1 0.1]);
axis([startTime endTime yLowerLim (a)])
set(gca, 'YTick', ticks, 'YTickLabel', value)
set(gca,'XTick',(startTime:1:endTime))
set(gca,'FontSize',20)
ylabel('Epicentral distance (km)')
xlabel('UTC times (hours)')

ax = gca;
yyaxis right
ax.YLim = ([yLowerLim a]);
set(ax, 'YTick', ticks, 'YTickLabel', staz_name)
box on

% SAVE THE TIME SERIES BOTH IN .FIG & .EPS FORMAT
if size(SatNum,2) == 1
    saveas(gcf,sprintf('TEC_timeseries_sat%d.fig',SatNum),'fig')
    print(gcf,'-depsc','-r0',sprintf('TEC_timeseries_sat%d',SatNum));
%     print(gcf,'-dpng','-r0','TEC_timeseries');
else 
    saveas(gcf,'TEC_timeseries_sat','fig')
    print(gcf,'-depsc','-r0','TEC_timeseries_sat');
end

%% Spectra analysis

%     nameDOY.PRNxx(:,01)  = Time
%     nameDOY.PRNxx(:,02)  = LAT_IPP
%     nameDOY.PRNxx(:,03) = LON_IPP
%     nameDOY.PRNxx(:,04) = ELE
%     nameDOY.PRNxx(:,05) = AZIM
%     nameDOY.PRNxx(:,06) = EFM
%     nameDOY.PRNxx(:,07) = STEC
%     nameDOY.PRNxx(:,08) = VTEC_FILT
%     nameDOY.PRNxx(:,09) = STEC_FILT
%     nameDOY.PRNxx(:,10) = DISTANCE 


for ii = 1:length(id_staz)
  
    GNSSname = staz(id_staz(ii),1:end-8);
    for jj = 1:length(SatNum)
        
        if (isfield(eval(staz(id_staz(ii),1:end-4)),sprintf('PRN%d',SatNum(jj)))) 

            TIME = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.PRN%d(:,1)',SatNum(jj))));
            filtTECU = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.PRN%d(:,9)',SatNum(jj))));
            distVect = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.PRN%d(:,10)',SatNum(jj))));
            
            % downsample the vectors
            n=(size(TIME,1))/2;
        %     odd_ind = 1:2:2*n;
            odd_ind = 1:size(TIME,1);

            fs = 1/((TIME(2)-TIME(1))*3600);
            nfft = 2^nextpow2(length(filtTECU(odd_ind)));
            wind= 80; %80 is good for 15 or 30 sec data, for 120 sec data better to use 10 (Number of samples)

            %find distance
            distVect =distVect(odd_ind);
            idd=find (TIME(odd_ind) >= EventTimeHours);
            idd= idd(1);
            dist = distVect(idd);
                if ~isnan(dist)

                    spectroFabio_2v(TIME(odd_ind),filtTECU(odd_ind),GNSSname,SatNum(jj),EventTimeHours,fs,wind,nfft,dist,startTime,endTime,Tmin)

                    saveas(gcf,sprintf('%sspectro_%s_sat%d',spectroDirectory,GNSSname,SatNum(jj)),'fig')
                    print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_sat%d.eps',spectroDirectory,GNSSname,SatNum(jj)));
                    print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_sat%d.png',spectroDirectory,GNSSname,SatNum(jj)));

                    close all
                end

        end
    
    end
end


%% ==========================END OF THE SCRIPT=============================