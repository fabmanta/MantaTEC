%% CALCULATE TEC INTENSITY INDEX (TECII) OF AN EVENT BASED ON ENTIRE WEEK
% Author :Fabio Manta 
% Last update: April-2020

% Data centers GNSS 
% https://www.unavco.org/data/gps-gnss/data-access-methods/dai2/app/dai2.html#
% ftp://cddis.gsfc.nasa.gov/gnss/data/daily/
% ftp://gssc.esa.int/gnss/data/daily/
% ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/
% ftp://ftp.earthobservatory.sg/SugarData
% http://garner.ucsd.edu/pub/rinex/

clear all;close all;clc;
addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull/'
cd '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/tec_processing';
warning off
% =========================================================================
% Set up Directories and Event info

fileName = 'list_EQ_onlyOcean_3v.xlsx';
TablesDirectory = ('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/CID_vs_Uplift/data_table/');
data = readtable(fullfile(TablesDirectory,fileName));
NameFields= (data.Properties.VariableNames)';
first_run = 'n';
dosave = 'n';


lowFreq = 0.001;%Hz for meteorite 0.0008; Hz for EQ 0.001; Hz for volcs 0.002
highFreq = 0.01;%Hz for meteorite 0.016; Hz for EQ end volcanoes 0.01; 
ele_cutoff = 10; % Elevation cut-off in degrees from the horizon
eventType = 'EQ';  %Volcano or EQ

if first_run == 'y'
    
    evntName = 'Mw7.8'; %sinabung or Mw7.0
    dateString = '20101025';
    defaultTimeInterval = 30; %aquisition period in seconds if known
    lat_epi = -3.48	;%N-S;
    lon_epi = 100.08;%W-E;
   
    Tmin=10; Tmax=20;% Time range to be analyze
    EventTime = datenum(2010,10,25,14,42,00); %UTC 
    [Year,Month,Day,Hours,Minutes,Seconds]= datevec(EventTime);
    Year = num2str(Year);
    EventTimeHours = Hours+Minutes/60+Seconds/3600;
    oneHourBefore = 0.5;% in UTC time
    twoHourAfter = 2.5;%
    fourHourAfter = 4;%
else
    for nn = 1:size(data.date,1)
    evntName(nn,:) = sprintf('Mw%0.1f',data.Mw(nn)); %sinabung or Mw7.0
    dateString(nn,:) = sprintf('%s',datestr(data.date(nn),'yyyymmdd'));%'20030925';
    defaultTimeInterval(nn) = data.acqTime_sec_(nn); %aquisition period in seconds if known
    lat_epi(nn) = data.Lat(nn);%N-S;
    lon_epi(nn) = data.Lon(nn);%W-E;
    EventTime(nn) = datenum(data.date(nn).Year,data.date(nn).Month,data.date(nn).Day,19,50,06); %UTC 
    Year(nn,:) = num2str(data.date(nn).Year);
    [~,~,~,Hours,Minutes,Seconds] = datevec(data.time_UT_(nn));
    EventTimeHours(nn) = Hours+Minutes/60+Seconds/3600;
    oneHourBefore(nn) = round(EventTimeHours(nn))-1.5;% in UTC time
    twoHourAfter(nn) = round(EventTimeHours(nn))+1.5;%
    fourHourAfter(nn) = round(EventTimeHours(nn))+4;%
%     if twoHourAfter(nn) >= 24
%         twoHourAfter(nn) = 24;
%     end
%     if fourHourAfter(nn) >= 24
%         fourHourAfter(nn) = 24;
%     end
    
    end
end

altitutePPI = 300; %ionosphere altitude in km
alt_ion = altitutePPI *1000; %ionosphere altitude in meters
utc2local = nan; % Number of hours + UTC
global staz;

%%

for kk = 1:size(lat_epi,2)
    
    SITE_COORD = []; %coordinates of the GNSS stations 
    staz = [];
    
    rinexFolder = sprintf('%s_week',dateString(kk,:));
    outputFolder = sprintf('output_%s_week',dateString(kk,:));
    plotFolder = sprintf('plots_%s_week',dateString(kk,:));

    rinexDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/1_data/%s_data/%s/',eventType,rinexFolder);
    outputDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/4_output/%s_output/%s/',eventType,outputFolder);
    plotsDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/5_plots/%s_plots/%s/',eventType,plotFolder);

    % http://www.gnsscalendar.com/index.html?year=2018
    % SP3 files are in GPS week
    sp3 = dir(sprintf('%s*.sp3',rinexDirectory));
    sp3 = struct2cell(sp3);
    sp3fileName = cell2mat(sp3(1,:)'); %ftp://igs.ensg.ign.fr/pub/igs/products/

    ionex = dir(sprintf('%s*.%si',rinexDirectory,Year(kk,end-1:end)));
    ionex = struct2cell(ionex);
    ionexFileName = cell2mat(ionex(1,:)'); %ftp://ftp.aiub.unibe.ch/CODE/

    mkdir(outputDirectory)
    mkdir(plotsDirectory)
    %% Calculate all parameters 
    % It returns a stroucture .mat file for each station with respect to all sat 

    if first_run == 'y'

        rinex211_processor(rinexDirectory,sp3fileName,ionexFileName,outputDirectory,...
                           Year,defaultTimeInterval,lon_epi,lat_epi,altitutePPI,alt_ion,Tmin,Tmax,SITE_COORD,highFreq,lowFreq,ele_cutoff);
    end

    %% Load data of the Day of the event
    cd(outputDirectory)
    load('staz_coordinates.mat')
    load_file_tec
    cd(plotsDirectory)
    %% calculate TECII and plot spectrogram
    % Plot of the spectrogram on the week of the event

    %     staz_NAME.PRNxx(02)  = LAT
    %     staz_NAME.PRNxx(03) = LON
    %     staz_NAME.PRNxx(04) = ELE
    %     staz_NAME.PRNxx(05) = AZIM
    %     staz_NAME.PRNxx(06) = EFM
    %     staz_NAME.PRNxx(07)  = TEC
    %     staz_NAME.PRNxx(08)  = TEC_V_FILT
    %     staz_NAME.PRNxx(09)  = TEC_FILT
    %     staz_NAME.PRNxx(10)  = DISTANCE

    n_satel = data.sat(kk); % 9, 15, 18, 21, 26, 27
    dt =1800; %time interval in X axis (seconds)



    for jj = 1:size(n_satel,2) %loop Satellite
        figure('position',[0,0,650,1000])
        set(gcf,'color','w');       

        for ii = 1:size(staz,1) %loop GNSS stations


            GNSSname = staz(ii,1:end-8);

            TIME = eval(strcat(staz(ii,1:end-4),sprintf('.PRN%d(:,1)',n_satel(jj))));
            filtTECU = eval(strcat(staz(ii,1:end-4),sprintf('.PRN%d(:,9)',n_satel(jj))));
            distVect = eval(strcat(staz(ii,1:end-4),sprintf('.PRN%d(:,10)',n_satel(jj))));
            filtTECU = signal_attenuation(filtTECU);
            elevVect = eval(strcat(staz(ii,1:end-4),sprintf('.PRN%d(:,4)',n_satel(jj))));
            filtTECU = vertical_tecu_calculator_2v(filtTECU,elevVect,alt_ion,n_satel(jj));    


            %downsample the vectors
            n=(size(TIME,1))/2;
        %     odd_ind = 1:2:2*n;
            odd_ind = 1:size(TIME,1);

            fs = 1/defaultTimeInterval(kk);
            nfft = 2^nextpow2(length(filtTECU(odd_ind)));

            if defaultTimeInterval(kk) == 30 || defaultTimeInterval(kk) == 15 || defaultTimeInterval(kk) == 20
                wind= 80; %80 is good for 15 or 30 sec data, for 120 sec data better to use 10 (Number of samples)
            elseif defaultTimeInterval(kk) == 120
                wind= 10; %80 is good for 15 or 30 sec data, for 120 sec data better to use 10 (Number of samples)
            end

            %find distance
            distVect =distVect(odd_ind);
            idd=find (TIME(odd_ind) >= EventTimeHours(kk));
            idd= idd(1);
            dist = distVect(idd);


            GNSS = upper(GNSSname);
            [filtTECU_ungapped] = fillmissing(filtTECU,'previous');
            filtTECU_ungapped(isnan(filtTECU_ungapped)) = 0;

            Hours = floor(EventTimeHours(kk));
            fract = (EventTimeHours(kk)-Hours);
            Minutes = floor(fract*60);
            Seconds = ((fract*60)- Minutes)*60;
            EventTimeSecs = Hours*3600+ Minutes*60 + Seconds;
            TminSec = TIME(1)*3600;

            %calculatre spectrogram
            %     tmpfiltTECU = fillmissing(tmpfiltTECU,'previous');
            [S,F,T,P] = spectrogram(filtTECU_ungapped,wind,wind-1,nfft,fs,'yaxis'); 

%             P2 = abs(P+eps)./400;
            pdb.spectro = 10*log10(P/1e-12); %convert power measurements in decibels
            pdb.F = F;
            pdb.T =T+TminSec;
            name= sprintf('pdb%d',ii);
            assignin('base', eval('name'), pdb);   

                    sub = subplot(7,1,ii);
                    box on
                %         set(sub,'Position',[0.13 0.85-(0.17*(ii-1)) 0.77 0.13]);
                    pcolor(pdb.T,pdb.F*1.e+3,pdb.spectro)
                    shading flat
%                     set(gca, 'CLim', [80 110]);%80 100 86 110
                    set(gca, 'CLim', [100 130]);%80 100 86 110
                    set(gca, 'YTick',[0 2 5 10]);
                    if defaultTimeInterval(kk) == 120 
                        ylim([0 5])
                    else
                        ylim([0 10])
                    end
                    xlim([oneHourBefore(kk)*3600, twoHourAfter(kk)*3600])

                    set(gca,'XTick',(oneHourBefore(kk)*3600:dt:twoHourAfter(kk)*3600))
                    set(gca,'XTickLabel',[])
                    set(gca, 'TickDir', 'out')
                    colormap(hot)
                    if ii==7
                        hold on
                        v =axis;
                        plot([EventTimeSecs EventTimeSecs],[0, 10],'w','LineWidth',1.5)
                        plot([EventTimeSecs+480 EventTimeSecs+480],[0, 10],'--w','LineWidth',1)
                %             plot([DeltaTimeSec+2400 DeltaTimeSec+2400],[0, 10],'--w','LineWidth',1)
                        plot([v(1) v(2)],[2, 2],'--w','LineWidth',1) %Brunt–Väisälä frequency
                    end

                    hold on

        end
        %     set(gca,'XTick',(oneHourBefore*3600:dt:twoHourAfter*3600))
            set(gca, 'XTicklabel',datestr(seconds(oneHourBefore(kk)*3600:dt:twoHourAfter(kk)*3600),'HH:MM'));
            c =colorbar('eastoutside');
            set(c,'position',[0.9414 0.3624 0.0184 0.5648],'FontSize', 16)                 
            ylabel('Frequency (mHz)','FontSize', 8,'FontWeight','Bold');

            if (dosave == 'y')
                saveas(gcf,sprintf('%sspectro_%s_sat%d_vtec',plotsDirectory,GNSSname,n_satel(jj)),'fig')
    %             print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_sat%d.eps',plotsDirectory,GNSSname,n_satel(jj)));
                print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_sat%d_vtec.png',plotsDirectory,GNSSname,n_satel(jj)));
            end

    end
    %% Version 0 Calculate the TECII
%     openfig(fullfile(plotsDirectory,sprintf('spectro_%s_sat%d.fig',GNSSname,n_satel)))
    T1 = 0;T2= 2;
    nhours = T1+T2;
    % nhours = 2;

    boundary1 = find(pdb1.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb1.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb1.spectro(find(pdb1.spectro<0))= nan;
    [tsh1]= max(max(pdb1.spectro(:,boundary1:boundary2)));
    subplot(7,1,1);plot([pdb1.T(boundary1) pdb1.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb1.T(boundary2) pdb1.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb2.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb2.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb2.spectro(find(pdb2.spectro<0))= nan;
    [tsh2]= max(max(pdb2.spectro(:,boundary1:boundary2)));
    subplot(7,1,2);plot([pdb2.T(boundary1) pdb2.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb2.T(boundary2) pdb2.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb3.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb3.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb3.spectro(find(pdb3.spectro<0))= nan;
    [tsh3]= max(max(pdb3.spectro(:,boundary1:boundary2)));
    subplot(7,1,3);plot([pdb3.T(boundary1) pdb3.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb3.T(boundary2) pdb3.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb4.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb4.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb4.spectro(find(pdb4.spectro<0))= nan;
    [tsh4]= max(max(pdb4.spectro(:,boundary1:boundary2)));
    subplot(7,1,4);plot([pdb4.T(boundary1) pdb4.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb4.T(boundary2) pdb4.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb5.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb5.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb5.spectro(find(pdb5.spectro<0))= nan;
    [tsh5]= max(max(pdb5.spectro(:,boundary1:boundary2)));
    subplot(7,1,5);plot([pdb5.T(boundary1) pdb5.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb5.T(boundary2) pdb5.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb6.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb6.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb6.spectro(find(pdb6.spectro<0))= nan;
    [tsh6]= max(max(pdb6.spectro(:,boundary1:boundary2)));
    subplot(7,1,6);plot([pdb6.T(boundary1) pdb6.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb6.T(boundary2) pdb6.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb7.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb7.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb7.spectro(find(pdb7.spectro<0))= nan;
    [tsh7]= max(max(pdb7.spectro(:,boundary1:boundary2)));
    subplot(7,1,7);plot([pdb7.T(boundary1) pdb7.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb7.T(boundary2) pdb7.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    MBL = (tsh1+tsh2+tsh3+tsh4+tsh5+tsh6)/6; 
    PSD_max = tsh7;

    ITPI_v0 = ((PSD_max/MBL)-1)*100;  % TEC Intesity Index in percentage


    annotation('textbox',[0.15,0.85,0.1,0.1],'String',sprintf('ITPI V0=%0.2f%s',ITPI_v0,'%'),'FontSize',16,'EdgeColor','none','FontWeight','Bold','Color','Green')
    %set color bar
%     set(gca, 'CLim', [80 110]);
%     oldColorMap = colormap(hot);
%     unitbar=(97-80);
%     postick = ((MBL-80)*size(oldColorMap,1))/unitbar;
%     newColorMap = oldColorMap;
%     newColorMap(round(abs(postick)),:) = [0 1 0];
%     newColorMap(round(abs(postick)-1),:) = [0 1 0];
%     colormap(newColorMap);

%     if (dosave == 'y')
%         saveas(gcf,sprintf('%sspectro_%s_sat%d_ITPI_v0_%d',plotsDirectory,GNSSname,n_satel,nhours),'fig')
%     %     print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_sat%d_TECII_v0_%d.eps',plotsDirectory,GNSSname,n_satel,nhours));
%         print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_sat%d_ITPI_v0_%d.png',plotsDirectory,GNSSname,n_satel,nhours));
%     end

    % Version 1 Calculate the TECII
    % T1 = 0;T2= 3;
    % nhours = T1+T2;
    % 
%     openfig(fullfile(plotsDirectory,sprintf('spectro_%s_sat%d.fig',GNSSname,n_satel)))

    boundary1 = find(pdb1.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb1.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb1.spectro(find(pdb1.spectro<0))= nan;
    [tsh1]= nanmean(max(pdb1.spectro(:,boundary1:boundary2)));
    subplot(7,1,1);plot([pdb1.T(boundary1) pdb1.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb1.T(boundary2) pdb1.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb2.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb2.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb2.spectro(find(pdb2.spectro<0))= nan;
    [tsh2]= nanmean(max(pdb2.spectro(:,boundary1:boundary2)));
    subplot(7,1,2);plot([pdb2.T(boundary1) pdb2.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb2.T(boundary2) pdb2.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb3.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb3.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb3.spectro(find(pdb3.spectro<0))= nan;
    [tsh3]= nanmean(max(pdb3.spectro(:,boundary1:boundary2)));
    subplot(7,1,3);plot([pdb3.T(boundary1) pdb3.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb3.T(boundary2) pdb3.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb4.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb4.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb4.spectro(find(pdb4.spectro<0))= nan;
    [tsh4]= nanmean(max(pdb4.spectro(:,boundary1:boundary2)));
    subplot(7,1,4);plot([pdb4.T(boundary1) pdb4.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb4.T(boundary2) pdb4.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb5.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb5.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb5.spectro(find(pdb5.spectro<0))= nan;
    [tsh5]= nanmean(max(pdb5.spectro(:,boundary1:boundary2)));
    subplot(7,1,5); plot([pdb5.T(boundary1) pdb5.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb5.T(boundary2) pdb5.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb6.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb6.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb6.spectro(find(pdb6.spectro<0))= nan;
    [tsh6]= nanmean(max(pdb6.spectro(:,boundary1:boundary2)));
     subplot(7,1,6); plot([pdb6.T(boundary1) pdb6.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                plot([pdb6.T(boundary2) pdb6.T(boundary2)],[0, 10],'y','LineWidth',1.5)

    boundary1 = find(pdb7.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb7.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb7.spectro(find(pdb7.spectro<0))= nan;
    [tsh7]= max(max(pdb7.spectro(:,boundary1:boundary2)));
    subplot(7,1,7);plot([pdb7.T(boundary1) pdb7.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb7.T(boundary2) pdb7.T(boundary2)],[0, 10],'y','LineWidth',1.5)


    MBL = (tsh1+tsh2+tsh3+tsh4+tsh5+tsh6)/6; 
    PSD_max = tsh7;


    ITPI_v1 = ((PSD_max/MBL)-1)*100;  % in percentage


    annotation('textbox',[0.4,0.85,0.1,0.1],'String',sprintf('ITPI V1=%0.2f%s',ITPI_v1,'%'),'FontSize',16,'EdgeColor','none','FontWeight','Bold','Color','Blue')
    %set color bar
%     set(gca, 'CLim', [80 110]);
%     oldColorMap = colormap(hot);
%     unitbar=(110-80);
%     postick = ((MBL-80)*size(oldColorMap,1))/unitbar;
%     newColorMap = oldColorMap;
%     newColorMap(round(abs(postick)),:) = [0 1 0];
%     newColorMap(round(abs(postick)-1),:) = [0 1 0];
%     colormap(newColorMap);

%     if (dosave == 'y')
%     saveas(gcf,sprintf('%sspectro_%s_sat%d_ITPI_v1_%d',plotsDirectory,GNSSname,n_satel,nhours),'fig')
%     % print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_sat%d_TECII_v1_%d.eps',plotsDirectory,GNSSname,n_satel,nhours));
%     print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_sat%d_ITPI_v1_%d.png',plotsDirectory,GNSSname,n_satel,nhours));
%     end

    % Version 2 Calculate the TECII
%     openfig(fullfile(plotsDirectory,sprintf('spectro_%s_sat%d.fig',GNSSname,n_satel)))
    % T1 = 0;T2= 2;
    % nhours = T1+T2;

    boundary1 = find(pdb1.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb1.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb1.spectro(find(pdb1.spectro<0))= nan;
    [tsh1]= nanmedian(max(pdb1.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb2.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb2.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb2.spectro(find(pdb2.spectro<0))= nan;
    [tsh2]= nanmedian(max(pdb2.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb3.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb3.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb3.spectro(find(pdb3.spectro<0))= nan;
    [tsh3]= nanmedian(max(pdb3.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb4.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb4.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb4.spectro(find(pdb4.spectro<0))= nan;
    [tsh4]= nanmedian(max(pdb4.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb5.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb5.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb5.spectro(find(pdb5.spectro<0))= nan;
    [tsh5]= nanmedian(max(pdb5.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb6.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb6.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb6.spectro(find(pdb6.spectro<0))= nan;
    [tsh6]= nanmedian(max(pdb6.spectro(:,boundary1:boundary2)));

    boundary1 = find(pdb7.T >= EventTimeSecs-(3600*T1));boundary1 = boundary1(1);
    boundary2 = find(pdb7.T >= EventTimeSecs+(3600*T2));boundary2 = boundary2(1);
    pdb7.spectro(find(pdb7.spectro<0))= nan;
    [tsh7]= max(max(pdb7.spectro(:,boundary1:boundary2)));
    subplot(7,1,7);plot([pdb7.T(boundary1) pdb7.T(boundary1)],[0, 10],'y','LineWidth',1.5)
                   plot([pdb7.T(boundary2) pdb7.T(boundary2)],[0, 10],'y','LineWidth',1.5)


    MBL = (tsh1+tsh2+tsh3+tsh4+tsh5+tsh6)/6; 
    PSD_max = tsh7;

    ITPI_v2 = ((PSD_max/MBL)-1)*100;  % in percentage


    annotation('textbox',[0.7,0.85,0.1,0.1],'String',sprintf('ITPI V2=%0.2f%s',ITPI_v2,'%'),'FontSize',16,'EdgeColor','none','FontWeight','Bold','Color','Red')
    %set color bar
%     set(gca, 'CLim', [80 110]);
%     oldColorMap = colormap(hot);
%     unitbar=(110-80);
%     postick = ((MBL-80)*size(oldColorMap,1))/unitbar;
%     newColorMap = oldColorMap;
%     newColorMap(round(abs(postick)),:) = [0 1 0];
%     newColorMap(round(abs(postick)-1),:) = [0 1 0];
%     colormap(newColorMap);

    fprintf('ITPI_v0 = %f, ITPI_v1 = %f, ITPI_v2 = %f\n',ITPI_v0,ITPI_v1,ITPI_v2)

    if (dosave == 'y')
    saveas(gcf,sprintf('%sspectro_%s_sat%d_ITPI_vtec',plotsDirectory,GNSSname,n_satel),'fig')
%     print(gcf,'-depsc','-r0',sprintf('%sspectro_%s_sat%d_TECII_v2_%d.eps',plotsDirectory,GNSSname,n_satel,nhours));
    print(gcf,'-dpng','-r0',sprintf('%sspectro_%s_sat%d_ITPI_vtec.png',plotsDirectory,GNSSname,n_satel));
    end
%     close all
end


    %%
    
%     % compare elevation angles
%     nstaz = 'cetr3380';
%     sa= 10;
%     t= sprintf('staz_%s(%d).timeFixedFLAT',nstaz,sa);
%     e= sprintf('staz_%s(%d).elevation',nstaz,sa);
%     d= sprintf('staz_%s(%d).distance',nstaz,sa);
% 
%     figure;
%     subplot(2,1,1)
%     plot(eval(t),eval(e),'b')
%     hold on;
%     % plot(staz_luzz3370(23).timeFixedFLAT,staz_luzz3370(23).elevation,'g')
%     plot([EventTime, EventTime],[0, 90],'--','color',[.7 .7 .7]);
%     set(gca,'XTick',(oneHourBefore:datenum(0,0,0,0,30,0):fourHourAfter))
%     axis([oneHourBefore fourHourAfter+datenum(0,0,0,0,00,0) 0 90])
%     xlabel('Time (Local)','FontSize', 12,'FontWeight','Bold');
%     ylabel('Elevation angle (Km)','FontSize', 12,'FontWeight','Bold');
%     datetick('x','HH:MM','keeplimits', 'keepticks')
%     % legend({'XMIS-PRN2', 'XMIS-PRN15','Event time'});
% 
%     subplot(2,1,2)
%     plot(eval(t),eval(d),'b')
%     hold on;
%     % plot(staz_luzz3370(23).timeFixedFLAT,staz_luzz3370(23).distance,'g')
%     plot([EventTime, EventTime],[0, 800],'--','color',[.7 .7 .7]);
%     set(gca,'XTick',(oneHourBefore:datenum(0,0,0,0,30,0):fourHourAfter))
%     axis([oneHourBefore fourHourAfter+datenum(0,0,0,0,00,0) 0 800])
%     xlabel('Time (Local)','FontSize', 12,'FontWeight','Bold');
%     ylabel('Epicentral dist (Km)','FontSize', 12,'FontWeight','Bold');
%     datetick('x','HH:MM','keeplimits', 'keepticks')
%     % legend({'XMIS-PRN2', 'XMIS-PRN15','Event time'});
