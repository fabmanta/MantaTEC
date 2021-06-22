function [fid,epochflag,clockOffset,date,numSV,SVlist] = ...
rinexReadsObsBlockHead211(line,fid)
% Reads the metadata of a RINEX 2.11 observations block.
%
% ATTENTION: Ignores all data regarding events flagged with numbers
% greater than 1!!!
%
% Positioned in a RINEX 2.11 GNSS observations text file at the beginning
% of an observation record of type "EPOCH/SAT OR EVENT FLAG", reads its
% contents (the metadata - one or more lines - of a block of observations)
%
% Inputs:
% fid: Matlab identifier of an open RINEX 2.11 GNSS observations text file
%
% Outputs:
% success: 1 if function performs successfully, 0 otherwise
% epochflag: Rinex observations epoch flag, as follows:
% 0: OK
% 1: power failure between previous and current epoch
% From now on the "event flags":
% 2: start moving antenna
% 3: new site occupation
% 4: header information follows
% 5: external event (epoch is significant)
% clockOffset: value of the receiver clock offset. If not present in the
% metadata of the observations block (it's optional RINEX 2.10
% data)it is assumed to be zero. If not zero implies that epoch,
% code, and phase data have been corrected by applying
% realtime-derived receiver clock offset
% date: time stamp of the observations block. Six-elements column-vector
% as follows:
% year: four-digits year (eg: 1959)
% month: integers 1..12
% day: integers 1..31
% hour: integers 0..24
% minute: integers 0..60
% second: reals 0..60
% numSV: number of satellites for which observations were made and are
% stored in the observations block
% SVlist: cell column-vector of the observed satellites three-characters
% strings (system identifier + one or two-digits prn code); it has as
% many elements as the number of the satellites observed. The system
% identifiers are: 'G', 'R', 'S' and 'E' (note that blank system
% identifiers are stored in SVlist as 'G')
%
% António Pestana, November 2012 and March 2015
% Copyright (c) by António Pestana

eof = 0;
success = 1;
gotData = 0;
epochflag = [];
clockOffset = [];
date = [];
numSV = 0;
SVlist = {};
% line = fgetl(fid); % returns -1 if only reads EOF

if line == -1
    eof = 1;
    disp(['INFO (rinexReadsObsBlockHead211): End of observations '...
    'text file reached'])
    return
end

% The first thing to do: the reading of the epoch flag
epochflag = str2num(line(29));
gotData = gotData + 1;

while epochflag > 1 % case of an event flag
    linejump = str2num(line(30:32));
    disp(['WARNING (rinexReadsObsBlockHead211): observations event '...
    'flag encountered; ' num2str(linejump) ...
    ' lines were ignored.'])
    for count=1:linejump + 1
    line = fgetl(fid);
    end
    epochflag = str2num(line(29));
end

% Reads the time stamp of the observations block (6 numerical values)
date = cell2mat(textscan(line,'%f', 6));
year = date(1,1);
if (year > 79)&&(year <= 99)
    year = year + 1900;
elseif (year >= 0)&&(year <= 79)
    year = year + 2000;
else
    success = 0;
    disp(['ERROR (rinexReadsObsBlockHead211): observations block time '...
    'stamp format unrecognized (not RINEX 2.11)!'])
    return
end
date(1,1) = year;

% Gets the number of used satellites
numSV = str2num(line(30:32)); % number of used SV

% Gets the receiver clock offset. It's optional data!
clockOffset = 0;
if size(line,2) == 80
    clockOffset = str2num(line(69:80));
end

% Creats the list of the satellites used in the observations block
col = 33;
for sat=1:numSV
    if ismember(sat,13:12:999)
        line = fgetl(fid);
        col = 33;
    end
    aux = line(col:col+2);

    if aux(1) == ' ' % Case of blanks for GPS identifiers
        aux(1) = 'G';
    end
    
SVlist = [SVlist; aux];
col = col + 3;
end
%%%%%%%%%%%%%%%%%%%%%%%% end rinexReadsObsBlockHead211 %%%%%%%%%%%%%%%%%%%%%%%%
