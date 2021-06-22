function tgd=get_tgd_ionex(ionexfile)

% get_tgd_ionex     Get Transmiter Group Delay (TGD) from a ionexfile
%             tgd = get_tgd_ionex(ionexfile)
%
%       ionexfile = rinex nav file
%       tgd = n x 2 matrix (n = number of satellites in rinex file)
%             column 1 = PRN number
%             column 2 = tgd in seconds

tgd=[];
% message
disp(['Getting transmitter group delays from ' ionexfile '...']);

fid = fopen(ionexfile);
% read header
line = fgetl(fid);

while 1
 line = fgetl(fid); 
  if ~isstr(line), break, end
  if ~isempty(findstr(line,'PRN / BIAS / RMS')) && line(4)=='G'
         tmpline=line(5:27);
	 delay = sscanf(tmpline, '%f',[1 3]);
         tgd=[tgd;delay];
   end 

end
tgd(:,3)=[];
tgd(:,2)=tgd(:,2)*1e-9;
