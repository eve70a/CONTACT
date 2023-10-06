%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ temp_dif, temp_tol ] = parse_its(fname, idebug)

if (nargin<2 | isempty(idebug))
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

headers = strvcat(' Case ', 'Temp-dep:');

% parse the contents of the file

temp_dif = [NaN]; temp_tol = [NaN];

while(~feof(f))

   % skip lines until one of the headers is found (case-sensitive)

   [ s, iline, ihdr ] = find_header(f, iline, headers, 1, idebug);
   if (ihdr<0), break; end

   % Process the consecutive section

   if (~isempty(strfind(s, 'Case')))

      % Case     1

      icase = sscanf(s,' Case %d');
      if (idebug>=5), 
         disp(sprintf('...Case %d starts at line %d',icase,iline));
      end

      % initialize output-values for this case

      temp_dif(icase,:)=NaN;
      temp_tol(icase,:)=NaN;
   end

   if (~isempty(strfind(s, 'Temp-dep:')))

      % 1, Temp-dep: |Tk-Tk-1|, .001 |Tk| :        103.9   1.039E-02

      tmp = sscanf(s, '%d, Temp-dep: |Tk-Tk-1|, .001 |Tk| : %f %f');
      it = tmp(1);
      temp_dif(icase,it) = tmp(2);
      temp_tol(icase,it) = tmp(3);
   end

end
fclose(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Helper function for locating the next section

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Helper function for reading one line of input

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

