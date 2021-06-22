function ifb=get_ifb_fabio(LGU,tgd_tmp,emf_tmp)

% solves for interfrequency bias
% works for 2 sucessive LGUs


A = 40.3;           % ionosphere constant
f1 = 1.57542e9;     % L1 frequency
f2 = 1.2276e9;      % L2 frequency
c = 0.299792458e9;  % speed of light
l1 = c/f1;
l2 = c/f2;
%                                                                                                                  
% decim = 3;          % decimation factor
% data_min = 2;       % minimum number of data per epoch
                                                                                                                 
% multiplication factors
K = (A * (f1^2-f2^2)) / (f1^2*f2*c);

% correct LGU for tgd

LGU=LGU+f2*tgd_tmp;

M=[K./emf_tmp -f2*ones(length(emf_tmp),1)];
out=M\LGU;

testVar=(M*out-LGU);
if testVar(1)==0 & testVar(2)==0 & out(1)>0
 ifb=out(2);
% disp(out(1)/1e16) 
else 
 ifb=NaN;
end

