function [dist,disth]=distance_sat_station_calculator_3v(lon_epi,lat_epi,SIP_mat,altitutePPI,sv)

% lonep: longitude of the epicenter (degree E)
% latep: latitude of the epicenter (degree N)
% lon: longitude of point (degree E)
% lat: latitude of point (degree N)
% h: altitude of point (km)
%
% dist: epicentral distance at earth surface
% disth: epicentral distance at altitude h
% alpha: azimuth toward epicenter (plane surface hypothesis)

% conversion to radian epicentral coordinates
lon_epi=lon_epi*pi/180;
lat_epi=lat_epi*pi/180;
radius=6371; %in Km


for jj = 1:length(sv)

    eval(['lat = SIP_mat.PRN' num2str(sv(jj)) '(:,2);']);
    eval(['lon = SIP_mat.PRN' num2str(sv(jj)) '(:,3);']);
    lon=lon*pi/180;
    lat=lat*pi/180;

    xx=acos(sin(lat_epi).*sin(lat)+cos(lat_epi).*cos(lat).*cos((lon-lon_epi)));
    dist(:,jj) = radius.*xx;
    disth(:,jj) = (radius+altitutePPI).*xx;
end

