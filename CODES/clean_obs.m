function [L1n,L2n,C1n,P1n,P2n]=clean_obs(L1,L2,C1,P1,P2)

for n=1:size(P1,2)
	[P1n(:,n),wind_p1]=clean_phase(P1(:,n));
	[P2n(:,n),wind_p2]=clean_phase(P2(:,n));
	[C1n(:,n),wind_c1]=clean_phase(C1(:,n));
	[L1n(:,n),wind_l1]=clean_phase(L1(:,n));
	[L2n(:,n),wind_l2]=clean_phase(L2(:,n));
	ww=wind_p1.*wind_p2.*wind_c1.*wind_l2.*wind_l1;
	L1n(:,n)=L1n(:,n).*ww;
	L2n(:,n)=L2n(:,n).*ww;
	P1n(:,n)=P1n(:,n).*ww;
	P2n(:,n)=P2n(:,n).*ww;
	C1n(:,n)=C1n(:,n).*ww;

end % for n

