
function [ opts2 ] = plot2d(sol, opt2)

%
% function [ opts2 ] = plot2d([sol], [opt2])
%
% Make a plot of tangential tractions Px or Py along a slice with constant x
% or constant y.
%
%  sol   == structure with CONTACT results as returned by loadcase.
%
%  opt2, opts2 == (in/output) structures with options for plot2d.
%           A valid structure opts2 is returned when plot2d is called with
%           no input arguments.
%
%   - orie  = flag for orientation of selected slices, 
%               1 = grid column(s) with constant x, 2 = grid row(s), const y.
%                   Default: 2.
%   - pxory = direction of tractions to plot, 'x' or 'y' for Px resp. Py.
%             Default: 'x'.
%   - xslc  = x-coordinate(s) for a column-wise plot (orie=1). Note: the 
%             grid column(s) closest to xslc is/are used. Default: 0.0.
%   - yslc  = y-coordinate(s) for a row-wise plot (orie=2). Note: the grid
%             row(s) closest to yslc is/are used. Default: 0.0.
%   - negpn = flag for vertical range to be plotted.
%              1 = show positive traction bound  "fstat*pn"
%              0 = show positive and negative traction bounds
%             -1 = show negative traction bound "-fstat*pn"
%   - facpt = multiplication factor for negating Px/Py.
%              1 = show original Px/Py.
%             -1 = show -Px (pxory=='x') or -Py (pxory=='y').
%   - xlim  = range of x-coordinates plotted (orie=2). Default: [x0-dx:x1+dx].
%   - ylim  = range of y-coordinates plotted (orie=1). Default: [y0-dy:y1+dy].
%   - plim  = vertical range of plot. Defaults: negpn=1: [0-eps:pmax+eps],
%             negpn=0: [-pmax-eps:pmax+eps], negpn=-1: [-pmax-eps:0+eps].
%   - pn_linestyle = Matlab linestyle for traction bound (default '--')
%   - pt_linestyle = Matlab linestyle for tangential tractions (default '-o')
% 

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% construct a local struct "myopt" in which default values are filled in for
% all available options.

myopt = struct( ...
  'orie',    2, ...
  'pxory', 'x', ...
  'xslc',   0., ...
  'yslc',   0., ...
  'negpn',   1, ...
  'facpt',   1, ...
  'xlim',   [], ...
  'ylim',   [], ...
  'plim',   [], ...
  'pn_linestyle', '--', ...
  'pt_linestyle', '-o'  ...
);

% If the user has not supplied any arguments,
%    return default options as first and only output argument

if (nargin<1)
   opts2 = myopt;
   return
end

% If the user has not supplied an opt2-struct, use the default

if (nargin<2 | isempty(opt2))
   opt2 = myopt;
end

% Check whether user-supplied opt2-struct contains unknown options

useropts = fieldnames(opt2);
ierror = 0;
for i = 1:length(useropts)
   if (~isfield(myopt, useropts{i}))
      ierror = ierror + 1;
      if (isempty(sol) | nargout>=1)
         disp(sprintf('Removing unknown option "%s"',useropts{i}));
         opt2=rmfield(opt2,useropts{i});
      else
         disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
         if (ierror==1)
            disp(sprintf('You may use "opt2=plot2d([],opt2);" to remove the option',useropts{i}));
         end
      end
   end
end

% Overwrite all values in "myopt" with user-supplied values

myopts = fieldnames(myopt);
for i = 1:length(myopts)
   if (isfield(opt2, myopts{i}) & ~isempty(mygetfield(opt2,myopts{i})))
      myopt = setfield(myopt, myopts{i}, mygetfield(opt2,myopts{i}));
   end
end

% Return options myopt: completed with defaults, removed unknown fields

if (isempty(sol))
   opts2 = myopt;
   return
end

% Add defaults for missing components of sol - for convenience of C.library

if (~isfield(sol, 'x'))
   sol.x = sol.xl + ([1:sol.mx]-0.5)*sol.dx;
end
if (~isfield(sol, 'y'))
   sol.y = sol.yl + ([1:sol.my]-0.5)*sol.dy;
end
if (~isfield(sol,'fric')), sol.fric.fstat = 0.3; end
if (~isfield(sol, 'mu'))
   sol.mu = sol.fric.fstat * ones(size(sol.pn));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start actual processing to produce the requested plot
%

%------------------------------------------------------------------------
% determine indices islc corresponding to requested x or y-value(s)

islc = [];

if (myopt.orie == 1)          % plot tractions along grid column with const x

   if (~isfield(myopt,'xslc') | isempty(myopt.xslc))
      disp('ERROR: when plotting along grid columns, opt2.xslc may not be empty');
      return
   end
   for j = 1:length(myopt.xslc)
      [xmin,ix] = min(abs(sol.x-myopt.xslc(j)));
      if (abs(sol.x(ix)-myopt.xslc(j))>0.001)
         disp(sprintf('Using grid column ix=%d with x=%f, nearest to x=%f', ...
                                                 ix, sol.x(ix), myopt.xslc(j)));
      end
      islc(j) = ix;
   end

elseif (myopt.orie == 2)      % plot tractions along grid row with const y

   if (~isfield(myopt,'yslc') | isempty(myopt.yslc))
      disp('ERROR: when plotting along grid rows, opt2.yslc may not be empty');
      return
   end
   for j = 1:length(myopt.yslc)
      [ymin,iy] = min(abs(sol.y-myopt.yslc(j)));
      if (abs(sol.y(iy)-myopt.yslc(j))>0.001)
         disp(sprintf('Using grid row iy=%d with y=%f, nearest to y=%f', ...
                                                 iy, sol.y(iy), myopt.yslc(j)));
      end
      islc(j) = iy;
   end

else

   disp(['ERROR: invalid value for orientation: ',num2str(myopt.orie)]);
   return

end

%------------------------------------------------------------------------
% get information for horizontal axis for appropriate orientation

if (myopt.orie == 1)          % plot tractions along grid column with const x

%  determine range of indices ia0/1 with non-zero pn (px/py)

%  the following is not supported in Matlab 5.3:
%  ia0 = find( any(sol.pn(:,islc)>0, 2), 1, 'first')-1;
   ia0 = find( any(sol.pn(:,islc)>0, 2) ) - 1;
   ia1 = find( any(sol.pn(:,islc)>0, 2) ) + 1;
   if (isempty(ia0)), ia0=1; ia1=sol.my; end;
   ia0 = max(1, ia0(1)); 
   ia1 = min(sol.my, ia1(end));
   [ia0 ia1]

%  determine plot limits for horizontal axis

   if (~isfield(myopt,'ylim') | isempty(myopt.ylim))
      axlim = [ sol.yl-sol.dy, sol.yl+(sol.my+1)*sol.dy ];
   else
      axlim = myopt.ylim;
   end

   axcoor = sol.y;
   axname = 'Y-coordinate [mm]';

elseif (myopt.orie == 2)      % plot tractions along row y=const

%  determine range of indices ia with non-zero pn (px/py)

   ia0 = find( any(sol.pn(islc,:)>0, 1) ) - 1;
   ia1 = find( any(sol.pn(islc,:)>0, 1) ) + 1;
   if (isempty(ia0)), ia0=1; ia1=sol.mx; end;
   ia0 = max(1, ia0(1)); 
   ia1 = min(sol.mx, ia1(end));

%  determine plot limits for horizontal axis

   if (~isfield(myopt,'xlim') | isempty(myopt.xlim))
      axlim = [ sol.xl-sol.dx, sol.xl+(sol.mx+1)*sol.dx ];
   else
      axlim = myopt.xlim;
   end

   axcoor = sol.x;
   axname = 'X-coordinate [mm]';

else
   disp(['ERROR: invalid value for orientation: ',num2str(myopt.orie)]);
   return
end

%------------------------------------------------------------------------
% get tangential tractions of appropriate direction (px='x', py='y').

if (strcmp(myopt.pxory,'x'))          % plot tangential tractions px
   pt = sol.px;
   pname = 'P_x';
elseif (strcmp(myopt.pxory,'y'))      % plot tangential tractions py
   pt = sol.py;
   pname = 'P_y';
else
   disp(['ERROR: invalid value for p-direction: ',myopt.pxory]);
   return
end

% get selected slices
if (myopt.orie == 1)          % plot columns x=const, ix=islc
   trcbnd = sol.mu(ia0:ia1, islc) .* sol.pn(ia0:ia1, islc);
   pt     =     pt(ia0:ia1, islc);
else                        % plot rows y=const, iy=islc
   trcbnd = sol.mu(islc, ia0:ia1)' .* sol.pn(islc, ia0:ia1)';
   pt     =     pt(islc, ia0:ia1)';
end

%------------------------------------------------------------------------
% determine appropriate vertical-range plim for plot:

if (isfield(myopt,'plim') & ~isempty(myopt.plim))

%  vertical plot limits provided by user:

   pmin = myopt.plim(1);
   pmax = myopt.plim(2);

else

%  derive vertical plot limits from solution.

   pmax = max(max(trcbnd));

%  find nice round value for pstep such that pmax ~= 10*pstep
   pdec = 1e-9;
   pfac = [1, 2, 2.5, 5];
   found = 0;
   while (~found)
      for k = 1 : length(pfac)
         if (pmax <= 12*pfac(k)*pdec)
            found = 1;
            pstep = pfac(k)*pdec;
            break;
         end
      end
      if (~found)
         pdec = pdec * 10;
      end
   end

%  set limits according to expected sign of tangential tractions
   if (myopt.negpn>0)
      pmax = ceil(pmax / pstep)*pstep;
      pmin = -0.5*pstep;
   elseif (myopt.negpn==0)
      pmax = ceil(pmax / pstep)*pstep;
      pmin = -pmax;
   else
      pmin = -ceil(pmax / pstep)*pstep;
      pmax = 0.5*pstep;
   end
end

plim = [ pmin, pmax ];

%------------------------------------------------------------------------
% layout configuration parameters:

xtickheight   = 0.01*(plim(2)-plim(1));
xticklabelsep = 1*xtickheight;
ytickwidth    = 0.007*(axlim(2)-axlim(1));
yticklabelsep = 2*ytickwidth;

%------------------------------------------------------------------------
% plot normal traction, tangential traction

clf
hold on

% plot traction bound, positive and/or negated

if (myopt.negpn>=0)
   plot(axcoor(ia0:ia1),  trcbnd, myopt.pn_linestyle);
end
if (myopt.negpn<=0)
   plot(axcoor(ia0:ia1), -trcbnd, myopt.pn_linestyle);
end

% plot tangential tractions

plot(axcoor(ia0:ia1), myopt.facpt*pt, myopt.pt_linestyle);

% use scaling derived from solution or provided through arguments:

axis([axlim plim]);

% disable Matlab's axis

xtick = get(gca,'xtick');
ytick = get(gca,'ytick');
axis off

% plot own horizontal axis with ticks and tick-labels

plot(axlim,[0 0],'k-');
for x = xtick,
    if (abs(x)<1e-6), x=0; end
    plot([x x],[0 1]*xtickheight,'k');
    text(x, -xticklabelsep, num2str(x), 'verticalalignment','top', ...
                                        'horizontalalignment','center'); 
end

% plot vertical axis with ticks and tick-labels

plot([0 0],[pmin pmax]*1.1,'k-');
for y = ytick,
    if (abs(y)>1e-6)
       plot([-1 1]*ytickwidth,[y y],'k-');
       text(-yticklabelsep, y, num2str(y), 'horizontalalignment','right'); 
    end
end

% set labels on plot

text(axlim(2), -6*xticklabelsep, axname, ...
                 'verticalalignment','top', 'horizontalalignment','right'); 
text(yticklabelsep, 1.1*ytick(end)-0.1*ytick(end-1), ...
                              ['Tractions ',pname,', f_{stat}*P_n [N/mm^2]']);

if (nargout>=1)
   opts2 = myopt;
end
