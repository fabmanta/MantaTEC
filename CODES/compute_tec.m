function [iec] = compute_tec(L1,L2,C1,P1,P2,sv,tgd,ifb,emf_mat)

% compute_tec	Computes ionospheric electron content
%
%	L1, L2, C1, P1, P2 are observables matrices (output of read_rinexo)
%	tgd transmitter group delay from ephenerides, sv satellites available
%	iec = integrated electron content, unit = el/m^2
%	iec_dot = time derivative of integrated electron content
%       lgu = unbiased lg observable, unit = cycles
%
%               WARNING: IEC and LGU output here are corrected for the
%               phase ambiguity (using pseudorange data) and also for
%               the satellite (tgd) and receiver (ifb) interfrequency biases

% constant
A = 40.3;

% L1 and L2 frequencies
f1 = 1.57542e9;
f2 = 1.2276e9;
c = 0.299792458e9;
l1 = c/f1;
l2 = c/f2;

% multiplication factor
F = (f1^2*f2*c) / (A * (f1^2-f2^2) * 1e16);

% number of satellites
nsat = size(L1,2);

% number of observations
nobs = size(L1,1);

% find tgd for each satellite
for n=1:length(sv)
	zeil(n)=find(tgd(:,1)==sv(n));	
end % for n
tmp_tgd=tgd(zeil,2);

% compute PG and LG
if (isempty(P1))
  PG = (f2/c) * (P2-C1);
else
  PG = (f2/c) * (P2-P1);
end
LG = (L2 - (f2/f1)*L1);

% LG vary in opposite sense from PG
LG = -LG;

% initialize arrays
iec = zeros(nobs,nsat);
lgu = zeros(nobs,nsat);
ifb = [];

% for each satellite
for i = 1:nsat
   x = [0 ; PG(:,i)./(PG(:,i)+eps^2)];
   y = [PG(:,i)./(PG(:,i)+eps^2) ; 0];
   z = x - y;
   I = find(z);
   B = mean(nonzeros(PG(:,i)-LG(:,i)));
   
   % foreach data segment
   for j = 1:2:length(I)
      % get data segment
%       pg_tmp = PG(I(j):I(j+1)-1,i);
%       lg_tmp = LG(I(j):I(j+1)-1,i);
      % get rid of large values: few first or last values of pseudorange data
      % is often incorrect (very large). This should be improved or even done
      % earlier, for instance right after P1, P2, C1 are read.
%       J = find(abs(pg_tmp)<1000);
      if ~(isempty(I))
        % calculate phase ambiguity using pseudorange data
        % B should probably be calculated from high elevation data only
%         B =mean(pg_tmp(J)-lg_tmp(J));
        
        % calculate unbiased LG
        lgu(I(j):I(j+1)-1,i) = B + LG(I(j):I(j+1)-1,i);
        
        % calculate integrated electron content
%         iec(I(j):I(j+1)-1,i) = F * (+f2*tgd(i) +B+LG(I(j):I(j+1)-1,i) +f2*ifb);
        iec(I(j):I(j+1)-1,i) = F * (+f2*tmp_tgd(i) +B + LG(I(j):I(j+1)-1,i));

      else
        disp(['WARNING in get_ion: all pg_tmp data removed from this segment']);
      end
   end
   
   % calculate the receiver interfrequency bias
   for nn=1:2:size(lgu,1)-1
      ifbtmp = get_ifb_fabio(lgu(nn:nn+1,i),tmp_tgd(i),emf_mat(nn:nn+1,i));
      ifb = [ifb; ifbtmp];
   end
   
end

% correct data for receiver interfrequency bias
ifb_mean = mean(ifb(find(~isnan(ifb)))); %mean value for all sat
iec(find(iec)) = iec(find(iec))+ (F*f2*ifb_mean);


