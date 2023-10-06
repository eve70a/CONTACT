function [creep, force, subs_max] = parse_out3(fname, use_struct, idebug)
%
% function [ sol ] = parse_out3(fname, use_struct, [idebug])
%    or
% function [creep, force, subs_max] = parse_out3(fname, [use_struct], [idebug])
%
% Reads CONTACT output for basic contact cases (module 3) from out-file fname and 
% returns the overall values from it:
%   creep    = struct with creepages    [pen, cksi, ceta, cphi, ncon, nadh, nslip; nplast]
%   force    = struct with total forces [fn, fx, fy, mz, el.en, fricw, pmax, temp1, temp2]
%   subs_max = struct with maximum subsurface stresses and locations
%                                       [sighyd, sigvm, sigtr, sigma1, sigma3, sigxx, sigyy, sigzz]
% for all cases in the file.
%
% if use_struct=0, the outputs are returned in arrays instead of structs.

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(use_struct))
   use_struct=1; % return structs with values
end
if (nargin<3 | isempty(idebug))
   idebug=1; % display input when idebug>=2
end

% open file and read second line for checking that it is a .out-file.

if (isempty(strfind(fname,'.out')) & isempty(strfind(fname,'.ref_out')))
   fname=[deblank(fname),'.out'];
end
if (~exist(fname,'file'))
   disp(['ERROR: cannot find file ',fname])
   creep=[]; force=[];
   return
end
f=fopen(fname);

iline=0;
[s, iline] = read_line(f, iline, idebug);
[s, iline] = read_line(f, iline, idebug);
if (~strfind(s, 'detailed investigation of 3D'))
   disp('ERROR: this file doesn''t look like a CONTACT .out-file.');
   return;
end

% list the special lines that precede the output data

headers = strvcat('Case', ...
                  'KINEMATIC CONSTANTS', ...
                  'TOTAL FORCES, TORSIONAL', ...
                  'ELAST.EN.   FRIC.POWER', ...
                  'CONTACT STATISTICS', ...
                  'ABSMAX SIGHYD');

subs_keys = strvcat('SIGHYD', 'SIGVM', 'SIGTR', 'SIGMA1', 'SIGMA3', 'SIGXX', 'SIGYY', 'SIGZZ');
nsubs_keys = size(subs_keys,1);

% parse the contents of the file

pen = []; cksi = []; ceta = []; cphi = [];
fn = []; fx = []; fy = []; mz = []; elen = []; fric = []; pmax = []; temp1 = []; temp2 = [];
ncon = []; nadh = []; nslip = []; nplast = [];
subs_max = [];

while(~feof(f))

   % skip lines until one of the headers is found (case-sensitive)

   [ s, iline, ihdr ] = find_header(f, iline, headers, 1, idebug);
   if (ihdr<0), break; end

   % Process the consecutive section

   if (~isempty(strfind(s, 'Case')))

      % Case     1

      icase = sscanf(s,' Case %d');
      if (idebug>=5), disp(sprintf('...Case %d starts at line %d',icase,iline)); end

      % initialize output-values for this case

      pen(icase)=NaN;  cksi(icase)=NaN; ceta(icase)=NaN;  cphi(icase)=NaN;
      fn(icase)=NaN;   fx(icase)=NaN;   fy(icase)=NaN;    mz(icase)=NaN;
      elen(icase)=NaN; fric(icase)=NaN; pmax(icase)=NaN;  temp1(icase)=NaN;  temp2(icase)=NaN;
      ncon(icase)=NaN; nadh(icase)=NaN; nslip(icase)=NaN; nplast(icase)=NaN;
      subs_max(icase, 1:4, 1:nsubs_keys)=NaN;
   end

   if (~isempty(strfind(s, 'KINEMATIC CONSTANTS')))

      %  KINEMATIC CONSTANTS, SHIFT T=1, F=1:
      %      DT          VELOC     FX/FSTAT/FN   CETA        CPHI
      %      0.001       1.000     -0.6570       0.000       0.000    
      %  KINEMATIC CONSTANTS, ROLLING T=3, F=2:
      %       CHI          DQ        VELOC     FX/FSTAT/FN  FY/FSTAT/FN  CPHI
      %      0.000       1.000       1.000      0.2525E-01   0.000       0.000
      % get the creepages or forces

      if (idebug>=9), disp(sprintf('...Found kinematic constants at line %d',iline)); end
      [s1, iline] = read_line(f, iline, idebug);
      [s2, iline] = read_line(f, iline, idebug);
      if (strfind(s1,'CHI'))
         tmp = sscanf(s2, '%*f %*f %*f %f %f %f');
      else
         tmp = sscanf(s2, '%*f %*f %f %f %f');
      end
      if (strfind(s1, 'CKSI'))
         cksi(icase) = tmp(1);
      else
         fx(icase)   = tmp(1);
      end
      if (strfind(s1, 'CETA'))
         ceta(icase) = tmp(2);
      else
         fy(icase)   = tmp(2);
      end
      cphi(icase) = tmp(3);
   end

   if (~isempty(strfind(s, 'TOTAL FORCES, TORSIONAL')))
      if (idebug>=9), disp(sprintf('...Found total forces at line %d',iline)); end
      % just skip one line, next line should contain 'ELAST.EN.   FRIC.POWER')
      %                                           or 'ELAST.EN.   FRIC.WORK')
      [s , iline] = read_line(f, iline, idebug);
   end

   if (~isempty(strfind(s, 'ELAST.EN.  FRIC.')) | ~isempty(strfind(s, 'ELAST.EN.   FRIC.')))
      % T=1, F=2:
      % TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER
      %     FN          FX          FY          MZ         ELAST.EN.  FRIC.POWER
      %     25.00      0.1885       0.000       0.000       4.736       0.000
      %     FN/G       SHIFT X      SHIFT Y    APPROACH    MAX(T1)     MAX(T1)
      %     25.00      -2.708E-03  -1.106E-15  6.492E-03    637.0       637.0

      % get the creepages or forces

      if (idebug>=9), disp(sprintf('...Found contact patch forces at line %d',iline)); end
      [s1, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s1, '%f %f %f %f %f %f');
      fn(icase)   = tmp(1);
      mz(icase)   = tmp(4);
      elen(icase) = tmp(5);
      fric(icase) = tmp(6);

      [s2, iline] = read_line(f, iline, idebug);
      [s3, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s3, '%f %f %f %f %f %f');
      if (~isempty([strfind(s2, 'CREEP X'), strfind(s2, 'SHIFT X')]))
         cksi(icase)  = tmp(2);
      else
         fx(icase)    = tmp(2);
      end
      if (~isempty([strfind(s2, 'CREEP Y'), strfind(s2, 'SHIFT Y')]))
         ceta(icase)  = tmp(3);
      else
         fy(icase)    = tmp(3);
      end
      if (~isempty(strfind(s2, 'APPROACH')))
         pen(icase)   = tmp(4);
      end
      if (~isempty(strfind(s2, 'PMAX')))        % >= v22.1
         pmax(icase)  = tmp(5);
      end
      if (~isempty(strfind(s2, 'MAX(T1,T2)')))  % >= v22.1
         temp1(icase) = tmp(6);
      elseif (~isempty(strfind(s2, 'MAX(T1)'))) % <= v21.1
         temp1(icase) = tmp(5);
         temp2(icase) = tmp(6);
      end
   end

   if (~isempty(strfind(s, 'CONTACT STATISTICS')))

      % CONTACT STATISTICS
      % N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,
      % POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION
      %    NPOT   NCON   NADH  NSLIP  INORM  ITANG
      %   3731   1286    919    367      1      1
      % or
      %    NPOT   NCON   NADH  NSLIP  NPLAST INORM  ITANG
      %      30     26     14      8      4      1      1

      % get the number of elements in contact/adhesion/slip

      if (idebug>=5)
         disp(sprintf('...Found contact statistics at line %d',iline));
      end
      [s , iline] = read_line(f, iline, idebug); % N:
      [s , iline] = read_line(f, iline, idebug); % POT:
      if (isempty(strfind(s, 'CON:'))), % support older format
         [s , iline] = read_line(f, iline, idebug);
      end
      [s , iline] = read_line(f, iline, idebug); % NPOT
      [s1, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s1, '%f %f %f %f %f %f');
      ncon(icase)  = tmp(2);
      nadh(icase)  = tmp(3);
      nslip(icase) = tmp(4);
      if (~isempty(strfind(s, 'NPLAST')))
         nplast(icase) = tmp(5);
      end
   end

   if (~isempty(strfind(s, 'ABSMAX SIGHYD')))

      % ABSMAX SIGHYD =    -0.853 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % MAX     SIGVM =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
      % MAX     SIGTR =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
      % MAX    SIGMA1 =     0.008 AT (X,Y,Z) = (   0.000,   0.000,   1.250)
      % MIN    SIGMA3 =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGXX =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGYY =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGZZ =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
 
      % process all lines containing a maximum value

      while(~feof(f) & ~isempty(strfind(s, 'AT (X,Y,Z)')))

         % get the key number
         key = strtrim(s(8:14));
         ikey = 1;
         while(ikey<nsubs_keys & ~strcmp(key, deblank(subs_keys(ikey,:))))
             ikey = ikey + 1;
         end
         if (~strcmp(key, deblank(subs_keys(ikey,:))))
            disp(sprintf('Internal error: unknown key "%s".', key));
            break
         end

         % get the values
         tmp = sscanf(s(17:end),' %f AT (X,Y,Z) = ( %f, %f, %f)');

         % store in subs_max
         subs_max(icase, 1:4, ikey) = tmp;

         % read next line from file
         [ s, iline ] = read_line(f, iline, idebug);
      end
   end

end
fclose(f);

if (use_struct)
   sol = struct;
   sol.creep = struct('pen',pen, 'cksi',cksi, 'ceta',ceta, 'cphi',cphi, ...
                      'ncon',ncon, 'nadh',nadh, 'nslip',nslip, 'nplast',nplast);
   sol.force = struct('fn',fn, 'fx',fx, 'fy',fy, 'mz',mz, 'elen',elen, ...
                      'pmax',pmax, 'temp1',temp1, 'temp2',temp2);

   sol.subs = struct();
   if (~isempty(subs_keys))
      for ikey = 1 : nsubs_keys
         key = deblank(lower(subs_keys(ikey,:)));
         sol.subs = setfield(sol.subs, key, subs_max(:,:,ikey));
      end
   end
   creep = sol;
else
   creep = [pen; cksi; ceta; cphi; ncon; nadh; nslip; nplast]';
   force = [fn; fx; fy; mz; elen; fric; pmax; temp1; temp2]';
end

end % function parse_out3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_header: Helper function for locating the next section

function [ s, iline, ihdr ] = find_header(f, iline, headers, case_sens, idebug)

% function [ s, iline, ihdr ] = find_header(f, iline, headers, case_sens, idebug)
% case_sens = 1: case-sensitive matching; 0: case-insensitive

num_hdr = size(headers,1);
if (~case_sens)
   headers = lower(headers);
end
s       = [];
is_key  = 0;

while(~is_key & ~feof(f))

   % read next line from input

   [s, iline] = read_line(f, iline, idebug);
   if (iline==-39), idebug=100; end

   % convert to lower-case when doing case-insensitive matching

   if (case_sens)
      scmp = s;
   else
      scmp = lower(s);
   end

   % check if it contains one of the headers

   for iheader = 1: num_hdr
      if (~isempty(strfind(scmp, deblank(headers(iheader,:)))))
         is_key = iheader;
         break
      end
   end
end

s = strtrim( s );
if (feof(f))
   ihdr = -1;
else
   ihdr = is_key;
end

end % function find_header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_line: Helper function for reading one line of input

function [ s, iline ] = read_line(f, iline, idebug)

% function [ s, iline ] = read_line(f, iline, idebug)

s = fgets(f);
iline = iline + 1; 

if (idebug>=2),
   % debug: print line without trailing \n (?)
   if (length(s)>=2 & double(s(end-1))==13)
      disp(sprintf('line %4d: %s', iline, s(1:end-2) ))
   elseif (length(s)>=1 & double(s(end))==10)
      disp(sprintf('line %4d: %s', iline, s(1:end-1) ))
   else
      disp(sprintf('line %4d: %s', iline, s))
   end
end

end % function read_line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
