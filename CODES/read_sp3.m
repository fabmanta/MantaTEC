function [Xs,Ys,Zs,Ts] = read_sp3(sp3file,sv);

% READ_SP3	Reads sp3 file and returns position vectors
%		Xs, Ys, Zs (in kilometers) for satellite sv
%		and corresponding time in secs of the current day.
%		[Xs,Ys,Zs,Ts] = read_sp3(sp3file,sv);

fid=fopen(sp3file,'r');

Xs=[]; Ys=[]; Zs=[]; Ts=[];

while 1
  line = fgetl(fid);
  if ~isstr(line), break, end;

  if (line(1)=='*')
    [t] = sscanf(line, '%c %d %d %d %d %d %f');
  end

  if (line(1)=='P')
     tmp_sv=str2num(line(3:4));
     line=line(5:length(line));	
     [p] = sscanf(line, '%f %f %f %f', [1 4]);
     
     if (tmp_sv == sv)
       Xs = [Xs;p(1)]; Ys = [Ys;p(2)]; Zs = [Zs;p(3)]; Ts = [Ts;t(5)*3600+t(6)*60+t(7)];
     end
  end
end

fclose(fid);

