function [fid, rinexHeader, gnssType, markerName, antDelta, apcoords,...
numOfObsTypes, typesOfObs, tFirstObs, tLastObs, tInterval, ...
timeSystem, numHeaderLines, clockOffsetsON, leapSec, eof] = ...
rinexReadsObsFileHeader211(file)
% Extracts relevant data from the header of a RINEX GNSS observations file.
%
% Analyzes the header of a RINEX GNSS observation file and extracts
% relevant data.
%
% Limited tests have been done using RINEX 2.11 observation files.
%
%
% Input:
% file: RINEX observation file
%
% Outputs:
% success: 1 if the reading of the RINEX observations file seems to be
% successful, 0 otherwise
% warnings: 1 if the reading of the RINEX observations was done with
% warnings, 0 otherwise
% rHeader: cell column-vector containing the following data:
% rinexVersion: RINEX version number; string: '' if not specified
% rinexType: RINEX file type; char
% gnssType: GNSS system of the satellites observed; can be 'G', 'R',
% 'S', 'E' or 'M' that stand for GPS, GLONASS, Geostationary
% signal payload, GALILEO or Mixed (satellites from various of
% the previous systems); char
% rinexProgr: name of the software used to produce de RINEX GPS nav
% file; '' if not specified
% rinexDate: date/time of the RINEX file creation; '' if not
% specified
% markerName: name of the antenna marker; '' if not specified
% antDelta: column vector ot the three components of the distance from
% the marker to the antenna, in the following order - up, east and
% north; reals; null vector by default
% numOfObsTypes: number of different observation types stored in the
% RINEX file; THIS IS CRITICAL DATA!
% typesOfObs: cell column-vector containing the observation types; each
% observation type is a two-character string, the first one (a
% capital letter) is an observation code and the second one (a digit)
% is a frequency code. THIS IS CRITICAL DATA!
%
% According to RINEX 2.11 these codes are:
%
% C: Pseudorange GPS: C/A, L2C
% Glonass: C/A
% Galileo: All
% P: Pseudorange GPS and Glonass: P code
% L: Carrier phase
% D: Doppler frequency
% S: Raw signal strengths or SNR values
% as given by the receiver for the
% respective phase observations (see comments of
% function rinexReadsObsBlock211)
%
% Frequency code
% GPS Glonass Galileo SBAS
% 1: L1 G1 E2-L1-E1 L1
% 2: L2 G2 -- --
% 5: L5 -- E5a L5
% 6: -- -- E6 --
% 7: -- -- E5b --
% 8: -- -- E5a+b --
%
% Observations collected under Antispoofing
% are converted to "L2" or "P2" and flagged
% with bit 2 of loss of lock indicator (LLI);
% read comments of function rinexReadsObsBlock211
%
% tFirstObs: time stamp of the first observation record in the RINEX
% observations file; column vector of reals [YYYY; MM; DD; hh; mm;
% ss.sssssss]; THIS IS CRITICAL DATA
% tLastObs: time stamp of the last observation record in the RINEX
% observations file; column vector of reals [YYYY; MM; DD; hh; mm;
% ss.sssssss]. NaN by default. THIS IS RINEX 2.11 OPTIONAL DATA
% tInterval: observations interval; seconds. NaN by default. THIS IS
% RINEX 2.11 OPTIONAL DATA
% timeSystem: three-character code string of the time system used for
% expressing tfirstObs; can be GPS, GLO or GAL; THIS IS CRITICAL DATA
% numHeaderLines: total number of lines of the header
% clockOffsetsON: receiver clock offsets flag. O if no realtime-derived
% receiver clock offset was applied to epoch, code and phase data (in
% other words, if the file only has raw data), 1 otherwise.
% 0 by default. THIS IS RINEX 2.11 OPTIONAL DATA
% leapSec: number of leap seconds since 6-Jan-1980. UTC=GPST-leapSec.
% NaN by default. THIS IS RINEX 2.11 OPTIONAL DATA
% eof: end-of-file flag; 1 if end-of-file was reached, 0 otherwise
%
% Based in the work of Kai Borre
% António Pestana, March 2015
% Copyright (c) António Pestana

fid = fopen(file,'rt');
eof = 0;
success = 1;
warnings = 0;
got_info = 0;
numHeaderLines = 0;
antDelta = [0; 0; 0];

clockOffsetsON = 0;
numLinTypObs = 0;
timeSystem = '';
tFirstObs = [0; 0; 0; 0; 0; 0];
tLastObs = NaN;
tInterval = NaN;
typesOfObs = {};
leapSec = NaN;
numOfObsTypes = 0;
rinexHeader = {};

while 1 % Gobbling the header

    numHeaderLines = numHeaderLines + 1;
line = fgetl(fid); % returns -1 if only reads EOF
if line == -1
    eof = 1;
    if got_info > 5
        success = 1;
    else
        fprintf(['Some important data regarding the GNSS '...
        'observations may be missing...\nProceed with caution.'])
        warnings = 1;
    end
        break
end


answer = strfind(line,'END OF HEADER'); % [] if the string isn't found

if ~isempty(answer) % the end of the header was found
    if got_info > 5
        success = 1;
        break
    else
        fprintf(['Some important data regarding the GNSS '...
        'observations may be missing...\nProceed with caution.'])
        warnings = 1;
    end
    break
end

if numHeaderLines == 1
    rinexVersion = strtrim(line(1:9));
    rinexType = line(21);
    if rinexType ~= 'O'
        disp('ERROR: the file is not a RINEX observations data file!')
        success = 0;
        fclose(fid);
        return
    end
    gnssType = line(41); % reads the GNSS system type
    if ~ismember(gnssType, [' ' 'G' 'R' 'S' 'E' 'M'])
        disp(['ERROR: "' gnssType '"' ' is an unrecognized satellite '...
        'system type.'])
        success = 0;
        fclose(fid);
        return
    end
    if strcmp(gnssType,' ')
        gnssType = 'G';
    end
    got_info = got_info + 1;
end

answer = strfind(line,'PGM / RUN BY / DATE');
if ~isempty(answer)
rinexProgr = strtrim(line(1:20));
rinexDate = strtrim(line(41:60));
got_info = got_info + 1;
end

answer = strfind(line,'MARKER NAME');
if ~isempty(answer)
markerName = strtok(line);
got_info = got_info + 1;
end

answer = strfind(line,'ANTENNA: DELTA H/E/N');
if ~isempty(answer)
    for k = 1:3
    [number, line] = strtok(line); % finds the substring containing
    % the deltas of the antenna
    % relative to the marker
    antDelta (k,1) = str2num(number);
    end
    got_info = got_info + 1;
end

answer = strfind(line,'APPROX POSITION XYZ');
if ~isempty(answer)
    for k = 1:3
    [number, line] = strtok(line); % finds the substring containing
    % the deltas of the antenna
    % relative to the marker
    apcoords (k,1) = str2num(number);
    end
    got_info = got_info + 1;
end

answer = strfind(line,'# / TYPES OF OBSERV');
if ~isempty(answer)
    numLinTypObs = numLinTypObs + 1; % one more line of observations types
    line = strtrim(line(1:60)); % deletes '# / TYPES OF OBSERV'
    
    if numLinTypObs == 1
    [nObs, line] = strtok(line);
    numOfObsTypes = str2num(nObs);
    
        for k = 1:min(9,numOfObsTypes)
        [obsType, line] = strtok(line);
            if size(obsType,2) ~= 2 || ~ismember(obsType(1),['C' 'P'...
                'L' 'D' 'S']) || ~ismember(obsType(2),['1' '2' '3'...
                '4' '5' '6' '7' '8'])
                disp(['ERROR (rinexReadsObsHeader211): ' obsType...
                ' is a not a standard RINEX 2.11 observation type!'])
                success = 0;
                fclose(fid);
                return
            end
        typesOfObs = [typesOfObs; obsType];
        got_info = got_info + 1;
        end
        
    else
        
        for k = 1:min(9,numOfObsTypes-(numLinTypObs-1)*9)
        [obsType, line] = strtok(line);
            if size(obsType,2) ~= 2 || ~ismember(obsType(1),['C' 'P'...
                'L' 'D' 'S']) || ~ismember(obsType(2),['1' '2' '3'...
                '4' '5' '6' '7' '8'])
                disp(['ERROR (rinexReadsObsHeader211): ' obsType...
                ' is a not a standard RINEX 2.11 observation type!'])
                success = 0;
                B6
                fclose(fid);
                return
            end
        typesOfObs = [typesOfObs; obsType];
        end
    end
end


answer = strfind(line,'TIME OF FIRST OBS');
if ~isempty(answer)
line = strtrim(line(1:60)); % deletes 'TIME OF FIRST OBS'
    for k = 1:6
    [tok, line] = strtok(line); % finds the substrings containing
    % the components of the time of the
    % first observation (YYYY; MM; DD;
    % hh; mm; ss.sssssss) and specifies
    % the Time System used in the
    % observations file (GPST, GLOT or
    % GALT)
        switch k
        case 1
        yyyy = str2num(tok);
        case 2
        mm = str2num(tok);
        case 3
        dd = str2num(tok);
        case 4
        hh = str2num(tok);
        case 5
        mnt = str2num(tok);
        otherwise
        ss = str2num(tok);
        end
    end
    
tFirstObs = [yyyy; mm; dd; hh; mnt; ss];
got_info = got_info + 1;
aux = strtok(line);
    switch aux
    case 'GPS'
    timeSystem = 'GPST';
    case 'GLO'
    timeSystem = 'GLOT';
    case 'GAL'
    timeSystem = 'GALT';
        otherwise
        switch gnssType
        case 'G'
        timeSystem = 'GPST';
        case 'R'
        timeSystem = 'GLOT';
        case 'E'
        timeSystem = 'GALT';
        otherwise
        fprintf(['CRITICAL ERROR (rinexReadsObsHeader211):\n'...
        'The Time System of the RINEX observations file '...
        'isn''t correctly specified!\n'])
        success = 0;
        fclose(fid);
        return
        end
    end
end

answer = strfind(line,'TIME OF LAST OBS'); % This is an optional record
if ~isempty(answer)
    for k = 1:6
        [tok, line] = strtok(line); % finds the substrings containing
        % the components of the time of the
        % first observation (YYYY; MM; DD;
        % hh; mm; ss.sssssss)
        switch k
        case 1
        yyyy = str2num(tok);
        case 2
        mm = str2num(tok);
        case 3
        dd = str2num(tok);
        case 4
        hh = str2num(tok);
        case 5
        mnt = str2num(tok);
        otherwise
        ss = str2num(tok);
        end
    end
    tLastObs = [yyyy; mm; dd; hh; mnt; ss];
end

answer = strfind(line,'INTERVAL'); % This is an optional record
if ~isempty(answer)
    tInterval = str2num(strtok(line));
end

answer = strfind(line,'RCV CLOCK OFFS APPL'); % This is an optional record!
if ~isempty(answer)
    if (strtok(line)=='0')
        clockOffsetsON = 0;
    elseif (strtok(line)=='1')
        clockOffsetsON = 1;
    else
        success = 0;
        disp(['ERROR (rinexReadsObsHeader211): unrecognized '...
        'receiver clock offsets flag!'])
        fclose(fid);
        return
    end
end

answer = strfind(line,'LEAP SECONDS'); % This is an optional record
if ~isempty(answer)
    leapSec = str2num(strtok(line));
end

end

if numOfObsTypes == 0||sum(tFirstObs) == 0
    success = 0;
    fprintf(['CRITICAL ERROR (rinexReadsObsHeader211)!\nTake a look '...
    'at the RINEX observations file %s'], file);
    fclose(fid);
    return
end

rinexHeader = {rinexVersion; rinexType; gnssType; rinexProgr; rinexDate};
% fclose(fid);
%%%%%%%%% end rinexReadsObsFileHeader211 %%%%%%%%%
