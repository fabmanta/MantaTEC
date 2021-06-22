function E = get_emf_fabio(el,R,Hb,Ht,H,option)

% GET_MF    Calculates mapping function:
%             el = vector of satellite elevation angles (degrees)
%             R = mean earth radius (m)
%             Hb = height of bottom of ionosphere (m)
%             Ht = height of top of ionosphere (m)
%             H = height of maximum electron density (peak F2)
%             option = 1/2  --> 2 different ways of calculating mapping function
%             E = vector of mapping functions

% convert elevations to radians
elr = el .* (pi/180.0);


% calculate ionospheric "thickness"
Hion = Ht - Hb;

% calculate mapping function, option 1
if (option == 1)
  L = sqrt((R+Ht)^2 - R^2.*cos(elr).^2) - sqrt((R+Hb)^2 - R^2.*cos(elr).^2);
  E = Hion ./ L;
% calculate mapping function, option 2
elseif (option == 2)
  Q = R / (R+H);
  beta = asin( Q .* cos(elr) );
  E = cos(beta);
end
