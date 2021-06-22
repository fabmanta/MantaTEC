EventTime = datenum(2021,03,04,19,28,31); %UTC 
[Year,Month,Day,Hours,Minutes,Seconds]= datevec(EventTime);
EventTimeHours = Hours+Minutes/60+Seconds/3600;

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

staz_name ='ahti'; %hlfj HOFN
doy= '063';
n_satel = 5;
t = eval(sprintf('%s%s0.PRN%d(:,1)',staz_name,doy,(n_satel)));
lat = eval(sprintf('%s%s0.PRN%d(:,2)',staz_name,doy,(n_satel)));
lon = eval(sprintf('%s%s0.PRN%d(:,3)',staz_name,doy,(n_satel)));
tec = eval(sprintf('%s%s0.PRN%d(:,8)',staz_name,doy,(n_satel)));
azim = eval(sprintf('%s%s0.PRN%d(:,5)',staz_name,doy,(n_satel)));
ele = eval(sprintf('%s%s0.PRN%d(:,4)',staz_name,doy,(n_satel)));
dist = eval(sprintf('%s%s0.PRN%d(:,10)',staz_name,doy,(n_satel)));
% tec = signal_attenuation(tec);


figure;plot(t,tec,'g','LineWidth',3)
hold on
plot([EventTimeHours,EventTimeHours] , [-1,1],'LineWidth',3)
xlabel('Time(UTC)')
% xlim([0 5]);
ylim([-1.5 1.5]);
ylabel('VTEC')
set(gca,'FontSize',18)
% gname


dtec =diff(tec);
