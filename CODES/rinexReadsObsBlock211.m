function [Obs, LLI, SS] =...
rinexReadsObsBlock211(fid, numSV, numOfObsTypes)
% Reads all the observations from a RINEX observations block.
%
% Positioned at the beginning of the line immediately after the header of the
% observations block, reads all the observations in this block of a RINEX
% observations file. This function is to be used after using function
% rinexReadsObsFileHeader211.
%
% Inputs:
% fid: Matlab file identifier of a Rinex observations text file
% numSV: number of satellites observed at the selected epoch
% numOfObsTypes: number of different observation types stored in the
% Rinex observations text file (obtained from the header of the RINEX
% observations file using function rinexReadsObsFileHeader211)
%
% Outputs:
% Obs: array [numSV x numObs] of observations; reals. Uses 0 for missing
% observations
% LLI: array [numSV x numObs] of "loss of lock indicator" numeric codes
% (binary 000 - 111, decimal 0 - 7, see RINEX standard). According to
% http://facility.unavco.org/software/teqc/faqs.html and the RINEX
% standard the following decimal values have the meaning:
% 0 - binary 000: OK
% 1 - binary 001: loss-of-lock for this SV and this frequency (L1 or
% L2) (cycle slip possible, original meaning of the LLI flag)
% 2 - binary 010: switch wavelength factor to opposite of last
% WAVELENGTH FACT L1/2 record setting for this SV and this
% frequency (L1 or L2)
% 4 - binary 100: anti-spoofing (A/S) is on for this SV or some other
% signal degradation/noise enhancement is in effect
%
% I think that all the 3-bits possible combinations can be used. For
% example, binary 101 (decimal 5) should flag A/S and loss-of-lock;
% 111 (decimal 7) should flag A/S and loss-of-lock and switch
% wavelength factor
%
% IMPORTANT: this function assigns -999 to this indicator when it is
% blank in the file (RINEX 2.11 standard states that blank
% means "not known")
%
% SS: array [numSV x numObs] of RINEX signal strength "normalized"
% values. Integer values from 1 to 9 plus 0 meaning "not known". The
% RINEX 2.11 standard states:
%
% Signal strength projected into interval 1-9:
% 1: minimum possible signal strength
% 5: threshold for good S/N ratio
% 9: maximum possible signal strength
% 0 or blank: not known, don't care
%
% IMPORTANT: this function assigns -999 to this indicator when it is
% blank in the file (RINEX 2.11 standard states that blank
% means "don't care")
%
%
% Based in the work of Kai Borre
% António Pestana, March 2015
% Copyright (c) António Pestana
success = 1;
eof = 0;
Obs = zeros(numSV, numOfObsTypes);
LLI = zeros(numSV, numOfObsTypes);
SS = zeros(numSV, numOfObsTypes);

for sat = 1:numSV
    lineNum = 0;
    
    for obsNum = 1:numOfObsTypes
        if ismember(obsNum,(1:5:numOfObsTypes))
            line = fgetl(fid); % reads one line of text
            if line == -1
                eof = 1;
                disp(['ERROR (rinexReadsObsBlock211): the end of the '...
                'observations text file was reached unexpectedly'])
                success = 0;
%                 return
            end
        lineNum = lineNum + 1;
        end
        
        charPos = ((obsNum -(lineNum - 1)*5)-1)*16 + 1;
        
        if (charPos > size(line,2)) % if tere are less then 4 observation...
            newObs = 0;
            newLLI = 0;
            newSS = 0;
        else
            newObs = str2num(line(charPos:charPos+13)); % reads one observation
            if isempty(newObs)
            newObs = 0;
            end

            if (charPos+14>= size(line,2))
              newLLI = -999; 
              newSS = -999;
            else
                newLLI = str2num(line(charPos+14)); % loss of lock indicator
                if isempty(newLLI)
                newLLI = -999;
                end

                newSS = str2num(line (charPos+15)); % signal strength
                if isempty(newSS)
                newSS = -999;
                end
            end
            
        end
        
        % Stores the data
        Obs(sat,obsNum) = newObs;
        LLI(sat,obsNum) = newLLI;
        SS(sat,obsNum) = newSS;
        
    end
end
%%%%%%%%% end rinexReadsObsBlock211.m %%%%%%%%%