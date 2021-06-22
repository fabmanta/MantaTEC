%% MAIN SCRIPT TO CALCULATE TEC FROM RINEX FILES
% Author :Fabio Manta 
% Last update: March-2021

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
addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull/'
cd '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/tec_processing';


warning off
% =========================================================================
% Set up Directories and Event info
global staz; % list of GNSS stations
global sv;

defaultTimeInterval = 30; %aquisition period in seconds
lowFreq = 0.001;%Hz for meteorite 0.0008; Hz for EQ 0.001; Hz for volcs 0.002
highFreq = 0.010;%Hz for meteorite 0.016; Hz for EQ end volcanoes 0.01; 
ele_cutoff = 15; % Elevation cut-off in degrees from the horizon


eventType = 'Volcano';  %'Volcano' or 'EQ' or 'Meteorite'
evntName = 'StVincent'; %'sinabung' or 'Mw7.0'
dateString = '20210411';
lat_epi = 13.33333;%N-S;
lon_epi = -61.183333;%W-E;
altitutePPI = 300; %ionosphere altitude in km
alt_ion = altitutePPI *1000; %ionosphere altitude in meters
utc2local = nan; % Number of hours + UTC
Tmin=0; Tmax=24;% Time range to be analyze
SITE_COORD = []; %coordinates of the GNSS stations 
EventTime = datenum(2021,04,11,0,00,00); %UTC 
[Year,Month,Day,Hours,Minutes,Seconds]= datevec(EventTime);
Year = num2str(Year);
EventTimeHours = Hours+Minutes/60+Seconds/3600;
oneHourBefore = 0;% in UTC time
twoHourAfter = 12;%
fourHourAfter = 7.5;%

rinexFolder = sprintf('%s_%s_%dkm',dateString,evntName,altitutePPI);
outputFolder = sprintf('output_%s_%s_%dkm',dateString,evntName,altitutePPI);
plotFolder = sprintf('plots_%s_%s_%dkm',dateString,evntName,altitutePPI);

% rinexDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/1_data/%s_data/%s/',eventType,rinexFolder);
% outputDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/4_output/%s_output/%s/',eventType,outputFolder);
% plotsDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/5_plots/%s_plots/%s/',eventType,plotFolder);
rinexDirectory = sprintf('/Users/fabio/Dropbox/StVincentSoufriere/%s_data/input/',eventType);
outputDirectory = sprintf('/Users/fabio/Dropbox/StVincentSoufriere/%s_data/output/',eventType);
plotsDirectory = sprintf('/Users/fabio/Dropbox/StVincentSoufriere/%s_data/plots/',eventType);

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

% PRN = G or R
%     nameDOY.PRNxx(:,01)  = Time
%     nameDOY.PRNxx(:,02)  = LAT_IPP
%     nameDOY.PRNxx(:,03) = LON_IPP
%     nameDOY.PRNxx(:,04) = ELE
%     nameDOY.PRNxx(:,05) = AZIM
%     nameDOY.PRNxx(:,06) = rawSTEC
%     nameDOY.PRNxx(:,07)  = DiffSTEC
%     nameDOY.PRNxx(:,08)  = filtVTEC
%     nameDOY.PRNxx(:,09)  = filtSTEC
%     nameDOY.PRNxx(:,10)  = DISTANCE

rinex211_processor_2v(rinexDirectory,sp3fileName,ionexFileName,outputDirectory,...
                   Year,defaultTimeInterval,lon_epi,lat_epi,altitutePPI,alt_ion,Tmin,Tmax,SITE_COORD,highFreq,lowFreq,ele_cutoff);

%% Load data of the Day of the event
cd(outputDirectory)
load_file_tec
load('staz_coordinates.mat')
cd(plotsDirectory)
% =========================================================================
%% HODOCHRONES

TECU_delta = 0.07;
save_fig = 'y';

id_staz= [1:size(staz,1)];
% id_staz= [13 103];
% id_staz = val';
% id_staz= id_st;

% SatNum = 1:size(name_sat_srt,1); % all sats;
SatNum = [5 22 30];
type = 2; % 1 = diffTEC 2=vtec

% figure('units','normalized','outerposition',[0 0 1 1]);
% figure
for sat_loop = 1:length(SatNum) %loop satellites
    
    figure
%     figure('units','normalized','outerposition',[0 0 1 2])
    set(gcf,'color','w');
 
    for staz_loop = 1: size(id_staz,2) %loop stations
            
        if (isfield(eval(staz(id_staz(staz_loop),1:end-4)),sprintf('%s',name_sat_srt(SatNum(sat_loop),:))))

            tmpTIME = eval(strcat(staz(id_staz(staz_loop),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(sat_loop),:)))); % related with Lat Lon
            tmpDistance = eval(strcat(staz(id_staz(staz_loop),1:end-4),sprintf('.%s(:,10)',name_sat_srt(SatNum(sat_loop),:))));

            if (type == 1)
                tmpTecu = eval(strcat(staz(id_staz(staz_loop),1:end-4),sprintf('.%s(:,7)',name_sat_srt(SatNum(sat_loop),:))));    
                tmpTecu = diff(tmpTecu);
                tmpTIME = tmpTIME(1:length(tmpTecu));
                tmpDistance = tmpDistance(1:length(tmpTecu));

%             alpha = (radius.*cosd(tmpELE))./(radius+PPI);
%             obliquity = sqrt(1-(alpha.^2));
%             tmpTEC = tmpTEC.*obliquity;

            elseif (type == 2)
                tmpTecu = eval(strcat(staz(id_staz(staz_loop),1:end-4),sprintf('.%s(:,8)',name_sat_srt(SatNum(sat_loop),:))));    
%                 tmpTecu = signal_attenuation(tmpTecu);
            end
            
            SatTrajectory(staz_loop) = cplot(tmpTIME,tmpDistance,tmpTecu,'o','markerfacecolor','flat','MarkerSize',5); %'linewidth',5);%
            hold on 
       
        end      
%             pause()
    end
        plot([EventTimeHours, EventTimeHours],[0, 1500],'--','color',[.7 .7 .7]);
        plot([EventTimeHours, EventTimeHours+0.5556],[0, 2000],'--','color',[1 0 0]);   %1 km/s
        plot([EventTimeHours, EventTimeHours+2.2222],[0, 2000],'--','color',[1 0 0]);   %1 km/s

        set(gca, 'CLim', [-TECU_delta,TECU_delta]);
        set(gca,'FontSize',20)
        colorbar;
        
        if (name_sat_srt(SatNum(sat_loop),1) == 'G')
            colormap(jet)
        else
            colormap(hot)
        end
        
        set(gca,'XTick',(oneHourBefore:1:fourHourAfter))
        axis([oneHourBefore fourHourAfter 0 1500])
        xlabel('Time','FontSize', 20);
        ylabel('Epicentral Range (Km)','FontSize', 20);
        title(sprintf('VTEC-%s',name_sat_srt(SatNum(sat_loop),:)))
        
end 

% SAVE HODOCHRONES

        if (save_fig == 'y')
                for n = 1:length(SatNum)
                    namesatellite = name_sat_srt(SatNum(n),:);
                    saveas(figure(n),sprintf('%sAll_staz_%s',plotsDirectory,namesatellite),'fig')
%                     print(figure(n),'-depsc','-r0',sprintf('%sAll_staz_%s',plotsDirectory,namesatellite));
                    print(gcf,'-dpng','-r0',sprintf('%sAll_staz_%s',plotsDirectory,namesatellite));
                end
        close all
        end
        
%% FIND CLOSER STATIONS TO SOURCE

% conversion to radian epicentral coordinates
lon_epi_rad=lon_epi*pi/180;
lat_epi_rad=lat_epi*pi/180;
radius=6371; %in Km


for jj = 1:length(staz)

    lat_staz = SITE_COORD(jj,1);
    lon_staz = SITE_COORD(jj,2);
    lon_staz_rad=lon_staz*pi/180;
    lat_staz_rad=lat_staz*pi/180;

    xx=acos(sin(lat_epi_rad).*sin(lat_staz_rad)+cos(lat_epi_rad).*cos(lat_staz_rad).*cos((lon_staz_rad-lon_epi_rad)));
    dist_staz_epi(:,jj) = radius.*xx; %distance Epicenter - GNSS receiver
end

id_st = find(dist_staz_epi <= 1000);
id_lat =  find(SITE_COORD(:,1) <= lat_epi + 1);

[val]=intersect(id_st,id_lat);
% val is the id of the station which are both <1000km from the source & Souther then the source 
        
%% TRAJECTORY MAP
%Ionospheric Pierce Point trajectory

% idx = ismember(staz,'usp10630.21o','row');
% lidx = find(idx);

% id_staz= [1:size(staz,1)];
% id_staz= [13 103];
id_staz = val';
% id_staz= id_st;

% SatNum = [3 6 16];
% SatNum = (name_sat_srt)';

delta_deg = 13;

figure('units','normalized','outerposition',[0 0 1 1]);
% figure;
lat_ex=[lat_epi-delta_deg lat_epi+delta_deg-5];
lon_ex=[lon_epi-delta_deg-5 lon_epi+delta_deg-5];

m_proj('Mercator','longtitudes',lon_ex,'latitudes',lat_ex);

% m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
m_grid('linest','none','box','fancy','tickdir','in','fontsize',16);
%for higher resolution use gshhs_i, for full resolution use gshhs_f
m_gshhs_h('patch',[1.0000    1.0000    0.8750]);
% m_gshhs_i('color','k');
% m_coast('color',[0 .6 0]);

hold on
[xEvent, yEvent]= m_ll2xy(lon_epi,lat_epi);
plot(xEvent,yEvent,'pr','MarkerSize',20,'MarkerFaceColor', 'r');
[xstaz, ystaz]= m_ll2xy(SITE_COORD(:,2),SITE_COORD(:,1));
plot(xstaz(id_staz),ystaz(id_staz),'ok','MarkerSize',4,'MarkerFaceColor', 'k'); 
text(xstaz(id_staz)+0.001,ystaz(id_staz)+0.001,upper(staz(id_staz,1:4)),'FontSize', 10);

headWidth = 3;
headLength = 3;

marker = [ 'o' 's' 'd' 'v'];
colore = ['m' 'g' 'b' 'c'];
 set(gca,'FontSize',16)
 set(gcf,'color','w');
 box on  

  
 for jj= 1:size(id_staz,2) 

         
        for ii = 1:length(SatNum) 
            
          hold on    

%             plot all station tecu
            if (isfield(eval(staz(id_staz(jj),1:end-4)),sprintf('%s',name_sat_srt(SatNum(ii),:))))

            tmpTIME = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(ii),:)))); % related with Lat Lon
            tmpLAT = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,2)',name_sat_srt(SatNum(ii),:))));  
            tmpLON = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,3)',name_sat_srt(SatNum(ii),:))));
            tmpElevation = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,4)',name_sat_srt(SatNum(ii),:))));  
            [tmpLON, tmpLAT] = m_ll2xy(tmpLON,tmpLAT);

            startingPosition = find(tmpTIME <= EventTimeHours-0);
            finalPosition = find(tmpTIME <= EventTimeHours+3);
            eventPosition = find(tmpTIME <= EventTimeHours);

          
                    if ~isempty (startingPosition) 
                        if startingPosition(end)==length(tmpLON)                        
                             text(tmpLON(end),tmpLAT(end),sprintf('%s',name_sat_srt(SatNum(ii),:)));
                        else
                            startingPosition = startingPosition(end);
                            finalPosition = finalPosition(end);
                            eventPosition = eventPosition(end);
                            
                            plot(tmpLON(startingPosition:finalPosition),tmpLAT(startingPosition:finalPosition),'Color',[.7 .7 .7])
                            if (size(SatNum,2)<=4)
                                cc = [colore(ii),marker(ii)];
                                plot(tmpLON(eventPosition),tmpLAT(eventPosition),cc,'MarkerFaceColor',colore(ii),'MarkerEdgeColor','k','MarkerSize',9)
                                text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('%s',name_sat_srt(SatNum(ii),:)));
                            else
                              plot(tmpLON(eventPosition),tmpLAT(eventPosition),'ko','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','k','MarkerSize',9)
                              text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('%s',name_sat_srt(SatNum(ii),:)));
                            end
%                             text(tmpLON(finalPosition),tmpLAT(finalPosition),sprintf('G%d',SatNum(ii)+1),'FontSize',6);
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
                            text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('%s',name_sat_srt(SatNum(ii),:)));
                        else
                            cc = [colore(1),marker(1)];    
                            plot(tmpLON(eventPosition),tmpLAT(eventPosition),cc,'MarkerFaceColor',colore(1),'MarkerEdgeColor','k','MarkerSize',9)
                            text(tmpLON(eventPosition),tmpLAT(eventPosition),sprintf('%s',name_sat_srt(SatNum(ii),:)));                    
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
%         
    saveas(gcf,'map','fig')
    print(gcf,'-depsc','-r0','map');
    print(gcf,'-dpng','-r0','map');
                       
%% FILT TEC timeseries
% Only station paired with satellite N in order of distance form the epicenter at the time of the event

save_fig = 'y';

% idx = ismember(staz,'norf0630.21o','row');
% lidx = find(idx)

% SatNum = [30];
% SatNum = (1:52);

% id_selection= id_st;
% id_selection= (1:length(staz));
% id_selection= [49 73 76 154];
id_selection= val';

staz_new = staz(id_selection,:);
dist_mat = nan(size(staz_new,1),length(SatNum));
xLowerLim = -0.2;

type = 2; % 1 = diffTEC 2=vtec
radius = 6371;
PPI = 300;

for ii = 1:size(staz_new,1)
     dist_staz =[];

    for jj = 1:length(SatNum)
        
        
            if (isfield(eval(staz_new((ii),1:end-4)),sprintf('%s',name_sat_srt(SatNum(jj),:))))

                tmpTIME = eval(strcat(staz_new((ii),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(jj),:)))); 
                tmpDistance = eval(strcat(staz_new((ii),1:end-4),sprintf('.%s(:,10)',name_sat_srt(SatNum(jj),:)))); 

                eventPosition = find(tmpTIME <= EventTimeHours);
                if isempty(eventPosition)
                    eventPosition = 1;
                else 
                    eventPosition = eventPosition(end);
                end

                  dist_staz =[dist_staz, round(tmpDistance(eventPosition))];
            else
                  dist_staz =[dist_staz, nan];
            end
    end
    dist_mat(ii,:) = dist_staz; 
end

numberSatellites = size(SatNum,2);
tmpDIST = reshape(dist_mat.',1,[]);
[value id_sort] = sort(tmpDIST);
ind = ceil(id_sort/numberSatellites);

figure
a=0;
ticks = [];        
staz_name = [];

for jj = 1:7%length(ind)
    
    id_sat = id_sort(jj) -(ind(jj)-1)*numberSatellites;
    tmpTIME = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(id_sat),:)))); 
    tmpELE = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.%s(:,4)',name_sat_srt(SatNum(id_sat),:)))); 

    if (type == 1)
        tmpTEC = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.%s(:,7)',name_sat_srt(SatNum(id_sat),:)))); 
        tmpTEC = diff(tmpTEC);
        tmpTIME = tmpTIME(1:length(tmpTEC));
        tmpELE = tmpELE(1:length(tmpTEC));

        alpha = (radius.*cosd(tmpELE))./(radius+PPI);
        obliquity = sqrt(1-(alpha.^2));
        tmpTEC = tmpTEC.*obliquity;
        
     elseif (type == 2)
        tmpTEC = eval(strcat(staz_new(ind(jj),1:end-4),sprintf('.%s(:,8)',name_sat_srt(SatNum(id_sat),:)))); 
%         tmpTEC = signal_attenuation(tmpTEC);       
    end
        hold on;
        ticks = [ticks, a];
        pairs = [staz_new(ind(jj),1:4),sprintf('-%s',name_sat_srt(SatNum(id_sat),:))];
        staz_name = [staz_name; pairs];
        plot(tmpTIME,tmpTEC+a,'r')
        a = a + 0.2;

end

yyaxis left
ylimt = get(gca,'Ylim');
plot([EventTimeHours, EventTimeHours],[ylimt(1), ylimt(2)],'--','color',[0.1 0.1 0.1]);   %1 km/s
axis([oneHourBefore fourHourAfter xLowerLim (a)])
set(gca, 'YTick', ticks, 'YTickLabel', value)
set(gca,'XTick',(oneHourBefore:1:fourHourAfter))
set(gca,'FontSize',20)
ylabel('Epicentral distance (km)')
xlabel('UTC times (hours)')

ax = gca;
yyaxis right
ax.YLim = ([xLowerLim a]);
set(ax, 'YTick', ticks, 'YTickLabel', staz_name)
box on

    if (save_fig == 'y')
        if (type == 1)
            saveas(gcf,'diffVTEC_timeseries','fig')
            %     print(gcf,'-depsc','-r0',sprintf('diffVTEC_timeseries_sat%d',SatNum));
            print(gcf,'-dpng','-r0','diffVTEC_timeseries');
        else
            saveas(gcf,'VTEC_timeseries','fig')
            %     print(gcf,'-depsc','-r0',sprintf('VTEC_timeseries_sat%d',SatNum));
            print(gcf,'-dpng','-r0','VTEC_timeseries');
        end
    end
%% Spectra analysis

% PRN = G or R
%     nameDOY.PRNxx(:,01)  = Time
%     nameDOY.PRNxx(:,02)  = LAT_IPP
%     nameDOY.PRNxx(:,03) = LON_IPP
%     nameDOY.PRNxx(:,04) = ELE
%     nameDOY.PRNxx(:,05) = AZIM
%     nameDOY.PRNxx(:,06) = rawSTEC
%     nameDOY.PRNxx(:,07)  = DiffSTEC
%     nameDOY.PRNxx(:,08)  = filtVTEC
%     nameDOY.PRNxx(:,09)  = filtSTEC
%     nameDOY.PRNxx(:,10)  = DISTANCE

% idx = ismember(staz,'ttsf0990.21o','row');
% lidx = find(idx);

% SatNum = [1 3 7 9 11 17 19 23];

% id_staz= (1:length(staz));
% id_staz= [1:size(staz,1)];
% id_staz= lidx;


for ii = 1:length(id_staz)
%   try  
    GNSSname = staz(id_staz(ii),1:end-8);

    for jj = 1:length(SatNum)
        
        if (isfield(eval(staz(id_staz(ii),1:end-4)),sprintf('%s',name_sat_srt(SatNum(jj),:))))

            tmpTIME = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(jj),:)))); % related with Lat Lon
            tmpTEC = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.%s(:,8)',name_sat_srt(SatNum(jj),:)))); 
            tmpDIST = eval(strcat(staz(id_staz(ii),1:end-4),sprintf('.%s(:,10)',name_sat_srt(SatNum(jj),:)))); 

            %downsample the vectors
            n=(size(tmpTIME,1))/2;
        %     odd_ind = 1:2:2*n;
            odd_ind = 1:size(tmpTIME,1);

            fs = 1/((tmpTIME(2)-tmpTIME(1))*3600);
            nfft = 2^nextpow2(length(tmpTEC(odd_ind)));
            wind= 80; %80 is good for 15 or 30 sec data, for 120 sec data better to use 10 (Number of samples)

            %find distance
            tmpDIST =tmpDIST(odd_ind);
            idd=find (tmpTIME(odd_ind) >= EventTimeHours);
            idd= idd(1);
            dist = tmpDIST(idd);
            spectroFabio_2v(tmpTIME(odd_ind),tmpTEC(odd_ind),GNSSname,name_sat_srt(SatNum(jj),:),EventTimeHours,fs,wind,nfft,dist,oneHourBefore,fourHourAfter,Tmin)

            saveas(gcf,sprintf('%sspectro_%s_%s',spectroDirectory,GNSSname,name_sat_srt(SatNum(jj),:)),'fig')
            print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_%s.eps',spectroDirectory,GNSSname,name_sat_srt(SatNum(jj),:)));
            print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_%s.png',spectroDirectory,GNSSname,name_sat_srt(SatNum(jj),:)));

            close all
        end
    end
end
%% Plot IPP map movie (GIF)

proj_type = 'Mercator'; % 'lambert'
% id_staz = (1:length(staz));
% SatNum = 1:size(name_sat_srt,1); % all sats;
% SatNum = [15];
filename = 'movie_3.gif';
initial_time = 10; % in hours
dt= 30;
timeStep = (0:dt:3600*3);
delta_deg =10;
lat_ex=[lat_epi-delta_deg lat_epi+delta_deg];
lon_ex=[lon_epi-delta_deg lon_epi+delta_deg];
type = 2; % 1 = diffTEC 2=vtec
radius = 6371;
PPI = 300;

Plot_movie_4v(proj_type,lat_ex,lon_ex,lat_epi,lon_epi,id_staz,name_sat_srt,SatNum,initial_time,dt,timeStep,radius,PPI,type,filename,outputDirectory,plotsDirectory)
