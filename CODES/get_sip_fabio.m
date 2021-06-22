function [SIP emf_mat] = get_sip_fabio(SAT,sv,alt_ion,apcoords,Tmin,Tmax)

% GET_SIP	Calculate sub-ionospheric points + associated information
%		SAT = structure containing sat positions (created by int_sp3)
%		sv = vector of PRN numbers (created by readrinex)
%		alt_ion = ionospheric height
%		apcoords = ECEF coordinates of GPS site (XYZ, m) (created by readrinex)
%		Tmin = window start time (hours UT)
%		Tmax = window end time (hours UT)
%		lat_epi, lon_epi = location of event
%		dec, inc = declination and inclination of magnetic field (degres)
%
%		SIP = structure, for each field (= PRN):
%		  time in UT hour of current day
%		  latitude of SIPs in degrees
%		  longitude of SIPs in degrees
%		  elevation w.r.t GPS site in degrees
%		  azimuth w.r.t GPS site in degrees
%		  distance from event to SIP (m)
%		  azimuth from event to SIP (deg CW from N)
%		  angle between magnetic direction and sat-rec line-of-sight (at the rec)
%		
%		SIP = get_sip(SAT,sv,alt_ion,apcoords,Tmin,Tmax,lat_epi,lon_epi,dec,inc);
%

% message
disp(['Computing SIP info for alt_ion=' num2str(alt_ion/1000) ' km;']);

% define some constants (should be read from a file)
R = 6371000.0;      % mean Earth radius
Hb = 100000.0;      % height of bottom of ionosphere
Ht = 1200000.0;     % height of top of ionosphere
H = 300000.0;       % height of maximum electron density

% elevation mapping function option (see get_emf_fabio.m)
emf_opt = 2;

% get coordinates of GPS site
W=xyz2wgs([0 apcoords']);
lat_site = W(3)*pi/180;
lon_site = W(2)*pi/180;
alt_site = 0;

% initialize SIP structure
SIP = [];
emf_mat = [];

% compute SIP for each satellite
for i=1:length(sv)

   % read orbit data for satellite sv
   eval(['T  = SAT.PRN' num2str(sv(i)) '(1,:);']);
   eval(['Xs = SAT.PRN' num2str(sv(i)) '(2,:);']);
   eval(['Ys = SAT.PRN' num2str(sv(i)) '(3,:);']);
   eval(['Zs = SAT.PRN' num2str(sv(i)) '(4,:);']);

   if (~isempty(find(Xs)))
      % select data inside time window
      I = find(T>=Tmin & T<=Tmax);
      Xs = Xs(I);
      Ys = Ys(I);
      Zs = Zs(I);
      T = T(I)';  % don't forget this transpose...

      % matrix of satellite positions (ECEF, meters)
      S = [Xs' Ys' Zs'];
   
      % compute elevation angle, range, and azimuth w.r.t. site
      [azim,elev,hlen] = azelle(S,apcoords);


      % compute SIP coordinates
      alpha = cos(elev)./(1+alt_ion/R);
      alpha = acos(alpha) - elev;
      lat_sip = sin(lat_site).*cos(alpha) + cos(lat_site).*cos(azim).*sin(alpha);
      lat_sip = asin(lat_sip);
      lon_sip = (sin(alpha).*sin(azim)) ./ cos(lat_sip);
      lon_sip = asin(lon_sip) + lon_site;

%       % compute elevation mapping function
      emf = get_emf_fabio(elev.*180/pi,R,Hb,Ht,H,emf_opt);
      emf_mat =[emf_mat, emf];
      % fill out SIP structure
      sip_tmp = [T lat_sip*180/pi lon_sip*180/pi elev*180/pi azim*180/pi emf];
      field = ['PRN' num2str(sv(i))];
      SIP = setfield(SIP,field,sip_tmp);

   end
end

