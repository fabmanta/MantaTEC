% https://ds.iris.edu/gmap/#network=*&channel=BHZ&starttime=2010-04-14&endtime=2010-04-15&maxlat=67.4261&maxlon=-11.5938&minlat=61.1008&minlon=-26.8867&drawingmode=box&planet=earth

% LIST DATA CENTER: https://www.fdsn.org/datacenters/
% http://brunnur.vedur.is/
% http://orfeus-eu.org/webdc3/

% first send email request of mseed file to IRIS through BREQ_FAST service  
% https://ds.iris.edu/ds/nodes/dmc/manuals/breq_fast/
% http://ds.iris.edu/ds/nodes/dmc/forms/breqfast-request/
% http://service.iris.edu/fdsnws/station/docs/1/builder/


% decompres the file using ==> tar -zxvf Fabio_s_FIRST_Request.244207.tar.mseed


% ======================================================================================

FILENAME = 'Package_1611569917333.mseed';

addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull/readMSEED'
addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull'
TESTDATA = '/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/Rafael_stage/seismic/';

startTime = datenum(2010,04,14);
endTime = startTime + 23/24; % 1 day later 

filepath = fullfile(TESTDATA, FILENAME);
ds = datasource('miniseed', filepath );
chantag = ChannelTag('VI.ISKR..BHZ');

% w = waveform(filepath, 'seed');
w=waveform(ds, chantag, startTime, endTime);

fobj1 = filterobject('b',[0.5, 1], 2);
fobj2 = filterobject('b',[1, 2], 2);
fobj3 = filterobject('b',[2, 4], 2);
w5_1 = filtfilt(fobj1, w);
w5_2 = filtfilt(fobj2, w);
w5_3 = filtfilt(fobj3, w);
plot(w5_1,'xunit','date')
spectrogram(w5_1)

rsamobj1 = waveform2rsam(w5_1, 'median', 60.0*5);
rsamobj2 = waveform2rsam(w5_1, 'median', 60.0*5);
rsamobj3 = waveform2rsam(w5_3, 'median', 60.0*5);
rsamvector = [rsamobj1 rsamobj2 rsamobj3];

plot(rsamobj1);
legend('filt 0.5-1 Hz')
hold on;
% plot(rsamobj2);plot(rsamobj3);
% legend('1:default','2:default2','3:max-10s')
xlabel('Time(UTC)')
ylabel('RSAM(m/s)')
set(gca,'FontSize',18)

%%
FILENAME = 'BORG.II.mseed';

addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull/readMSEED'
addpath '/Users/fabio/Fabio_data/PostdocCNES_2019/0_MATLAB_codes/usefull'
TESTDATA = '/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/Rafael_stage/seismic/Fabio_s_FIRST_Request.244207';

startTime = datenum(2010,04,14);
endTime = startTime + 20/24; % 1 day later 


filepath = fullfile(TESTDATA, FILENAME);
ds = datasource('miniseed', filepath );
chantag = ChannelTag('II.BORG.00.BHZ');

w=waveform(ds, chantag, startTime, endTime);
% w = fillgaps(w,'interp');
% w = detrend(w);

fobj = filterobject('b',[0.5, 1], 2);
w5_1 = filtfilt(fobj, w);
% w5 = medfilt1(w5,100);

% M = movmedian(w5,3);
% plot(M,'xunit','date')

plot(w5_1,'xunit','date')
spectrogram(w5_1)


rsamobj1 = waveform2rsam(w5_1, 'median', 60.0*5)
% rsamvector = [rsamobj1];
plot(rsamobj1);
legend('filt 0.5-1 Hz')
hold on;
% rsamvector.plot();
% legend('1:default','2:default2','3:max-10s','4:median-600s')
xlabel('Time(UTC)')
ylabel('RSAM(m/s)')
set(gca,'FontSize',18)