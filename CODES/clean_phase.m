function [phase_clean, window]=clean_phase(phase)
% function phase_clean=clean_obs(phase)
%  cleans oscillations in the phase)
%  only segments longer than seg_leng survive
%	------> note: has to be expanded to all phase data as well
%
seg_len=120; % segment length in datapoints assuming 1/30s 0.5h =60
% max_disc=0.5e7;  %0.5e7; maximum difference between two phase values allowed
max_disc=0.5e6;  %0.5e7; maximum difference between two phase values allowed

	
A=[0;phase]; B=[phase;0];       % Calculate difference between following values
C=B-A; C=abs(C);  % phase(n)-phase(n-1) starts at n=2

jumps=(B~=0 & C<=max_disc );  % jumps =1 when C has a allowed value <max disc and the data is non zero

X=jumps(2:length(jumps)); Y=jumps(1:length(jumps)-1);
Z=X-Y;  % if Z=-1 jumpdown Z=1 jump up , assume no jump at first and last value
up=find(Z==1);	% from this point to the next, jumps goes from 0 -> 1
down=find(Z==-1);    % from this point to the next, jumps goes from 1 -> 0
up=up+1;

window=zeros(length(phase),1);
phase_clean=zeros(length(phase),1);

if ~isempty(up)
	if length(up)==length(down) & up(1)<=down(1)
		segments=[up down down-up+1]; % first one up second one down
	elseif  length(up)==length(down) & up(1)>down(1)
		segments=[down up up-down+1 ]; % first one down second one up
	elseif length(up)<length(down) & up(1)>=down(1)
		up=[1;up];
	        segments=[up down down-up+1];
	elseif length(up)>length(down)  & up(1)<down(1)
		down=down-1;
		down=[down;length(phase_clean)];	
		segments=[up down down-up+1];
	end %ifI

	segments=segments(find(segments(:,3)>=seg_len),:);
        
    for n=1:size(segments,1)
        phase_clean(segments(n,1):segments(n,2))=phase(segments(n,1):segments(n,2));
        window(segments(n,1):segments(n,2))=ones(length(segments(n,1):segments(n,2)),1);
    end % for n
end %if ~isempty