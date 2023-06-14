
function [ opt_out ] = plotstrs(sol, opt)

%
% function [ opt ] = plotstrs(sol, opt)
%
% Make a plot of quantities in the interiors of the contacting bodies.
%
%  sol   == structure with subsurface results as returned by loadstrs.
%
%  opt   == (input/output) structure with options for plotstrs.
%           A valid structure opt is returned when plotstrs is called with
%           no input arguments.
%
%   - block = selects a single block of subsurface points from the input
%   - field = quantity to plot,
%      'ux', 'uy', 'uz' for displacements,
%      'sighyd' or 'hydro' for the mean hydrostatic stress sighyd = I_1/3,
%      'sigvm' or 'mises' for the von Mises stress sigvm = sqrt(3*J_2),
%      'sigma1','sigma2','sigma3' for the principal stresses,
%      'sigtr' or 'tresca' for the maximum shear stress sigtr = sigma1-sigma3,
%      'sigxx', 'sigxy', 'sigxz', 'sigyy', 'sigyz', 'sigzz' for the
%               components of the stress tensor.
%     Default: 'mises'.
%   - dir   = orientation of selected slice, 'y' == plot 'Oxz'-plane (default)
%   - yslc  = y-coordinate(s) for an Oxz-plot (dir='y'). Default: 0.0.
%             Note: the grid column(s) closest to yslc is/are used.
%             Setting 'max' shows the maximum value over all slices.
%   - xslc, zslc : see yslc.
%   - addplot = clear/do not clear the figure before plotting;
%   - typplot = 'contour', 'contourf' (filled), 'surf'. Default: 'contourf'.
%   - cntrlvl = values at which contours are required. Default: 'auto'.
%   - clabel  = 'on' or [values]: add labels on contours
%   - scale   = 'linear' or 'log'. Note: log-scale takes absolute value of
%               the data, and works best with user-defined contour levels
%               (cntrlvl). Default: 'linear'.
%   - colormap = Matlab color-map to be used, e.g. 'cool', 'hot', 'jet', ...

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% determine version information on the Matlab graphics system

old_graphics = verLessThan('matlab', '8.4.0'); % 8.4 == R2014b
new_graphics = ~old_graphics;

% construct a local struct "myopt" in which default values are filled in for all available options.

myopt = struct( ...
   'field',    'mises', ...
   'dir',      'y', ...
   'xslc',     0., ...
   'yslc',     0., ...
   'zslc',     0., ...
   'addplot',  0, ... % clear plot (0) or add to exist.plot (1, experimental)
   'typplot',  'contourf', ...
   'cntrlvl',  'auto', ...
   'clabel',   'off', ...
   'scale',    'linear', ...
   'colormap', 'parula' ...
);
if (old_graphics)
   myopt.colormap = 'jet';
end

% If the user has not supplied any arguments,
%    return default options as first and only output argument

if (nargin<1)
   opt_out = myopt;
   return
end

% If the user has not supplied an opt-struct, use the default

if (nargin<2 | isempty(opt))
   opt = myopt;
end

% Check whether user-supplied opt-struct contains unknown options

useropts = fieldnames(opt);
ierror = 0;
for i = 1:length(useropts)
   if (~isfield(myopt, useropts{i}))
      ierror = ierror + 1;
      disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
      if (ierror==1)
         disp(sprintf('You may use "opt=rmfield(opt,''%s'');" to remove the option',useropts{i}));
      end
   end
end

% Overwrite all values in "myopt" with user-supplied values

myopts = fieldnames(myopt);
for i = 1:length(myopts)
   if (isfield(opt, myopts{i}) & ~isempty(mygetfield(opt,myopts{i})))
      myopt = setfield(myopt, myopts{i}, mygetfield(opt,myopts{i}));
   end
end

% Check whether requested field is valid

ok_fields=strvcat('ux', 'uy', 'uz', 'sighyd', 'hydro', 'sigvm', 'mises', ...
                  'sigma1', 'sigma2', 'sigma3', 'sigtr', 'tresca', ...
                  'sigxx', 'sigxy', 'sigxz', 'sigyy', 'sigyz', 'sigzz');
found=0;
for istr=1:size(ok_fields,1)
   if (strcmp(myopt.field, deblank(ok_fields(istr,:))))
      found=1;
   end
end
if (~found)
   disp(sprintf('Unknown field requested="%s"; Available options:',myopt.field));
   disp(ok_fields)
   return;
end

% Check whether clabel is valid

is_ok = (isempty(myopt.clabel) | strcmp(myopt.clabel,'auto') | strcmp(myopt.clabel,'on') | ... 
         strcmp(myopt.clabel,'off') | isnumeric(myopt.clabel) );
if (~is_ok)
   disp(sprintf('Unknown value clabel="%s"; Available options: on|off|auto or [values]',myopt.clabel));
   myopt.clabel = 'off';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start actual processing to produce the requested plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (myopt.addplot<=0)
   clf;
end
hold on;

if (strcmp(myopt.field,'ux'))

   show_2d_slice(sol, 'ux', myopt);
   title (['Displacement U_x in X-direction', get(get(gca,'title'),'string')]);
   hc=colorbar; log_colorbar(hc, opt);
   ylabel(hc, 'Displacement U_x [mm]');

end

if (strcmp(myopt.field,'uy'))

   show_2d_slice(sol, 'uy', myopt);
   title (['Displacement U_y in Y-direction', get(get(gca,'title'),'string')]);
   hc=colorbar; log_colorbar(hc, opt);
   ylabel(hc, 'Displacement U_y [mm]');

end

if (strcmp(myopt.field,'uz'))

   show_2d_slice(sol, 'uz', myopt);
   title (['Displacement U_z in Z-direction', get(get(gca,'title'),'string')]);
   hc=colorbar; log_colorbar(hc, opt);
   ylabel(hc, 'Displacement U_z [mm]');

end

if (strcmp(myopt.field,'sighyd') | strcmp(myopt.field,'hydro'))

   show_2d_slice(sol, 'sighyd', myopt);
   title (['Mean hydrostatic stress SIGHYD', get(get(gca,'title'),'string')]);
   hc=colorbar; log_colorbar(hc, opt);
   ylabel(hc, 'stress invariant \sigma_{hyd} [N/mm^2]');

end

if (strcmp(myopt.field,'sigvm') | strcmp(myopt.field,'mises'))

   show_2d_slice(sol, 'sigvm', myopt);
   title (['Von Mises stress \sigma_{vm}', get(get(gca,'title'),'string')]);
   hc=colorbar; log_colorbar(hc, opt);
   ylabel(hc, 'Von Mises stress \sigma_{vm} [N/mm^2]');

end

if (strcmp(myopt.field,'sigxx') | strcmp(myopt.field,'sigxy') | strcmp(myopt.field,'sigxz') | ...
    strcmp(myopt.field,'sigyy') | strcmp(myopt.field,'sigyz') | strcmp(myopt.field,'sigzz') )

   if (~isfield(sol, myopt.field))
      disp(['ERROR: no data available for ',myopt.field,'.']);
   else
      show_2d_slice(sol, myopt.field, myopt);
      title (['Stress component ',myopt.field,' ',get(get(gca,'title'),'string')]);
      hc=colorbar; log_colorbar(hc, opt);
      ylabel(hc, ['Stress component ',myopt.field,' [N/mm^2]']);
   end
end

if (strcmp(myopt.field,'sigma1') | strcmp(myopt.field,'sigma2') | strcmp(myopt.field,'sigma3'))

   if (~isfield(sol, myopt.field))
      disp(['ERROR: no data available for ',myopt.field,'.']);
   else
      show_2d_slice(sol, myopt.field, myopt);
      title (['Principal stress \sigma_',myopt.field(6),' ',get(get(gca,'title'),'string')]);
      hc=colorbar; log_colorbar(hc, opt);
      ylabel(hc, ['principal stress \sigma_',myopt.field(6),' [N/mm^2]']);
   end
end

if (strcmp(myopt.field,'sigtr') | strcmp(myopt.field,'tresca'))

   if (~isfield(sol, 'sigtr'))
      disp(['ERROR: no data available for ',myopt.field,'.']);
   else
      show_2d_slice(sol, 'sigtr', myopt);
      title (['Maximum shear stress \sigma_{tresca}',get(get(gca,'title'),'string')]);
      hc=colorbar; log_colorbar(hc, opt);
      ylabel(hc, ['maximum shear stress \sigma_{tresca} [N/mm^2]']);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_2d_slice(sol, field, opt)
% 
% function [ ] = show_2d_slice(sol, field, opt)
%
% generic routine for a 2d contour plot of 3d data
%   sol     - structure with subsurface results as returned by loadstrs
%   field   - the actual array to be plotted, or a field name within sol
%   opt     - structure with options as defined by plotstrs.
%

   % get field to be plotted

   nx  = sol.nx;
   ny  = sol.ny;
   nz  = sol.nz;

   if (ischar(field))
      eval(['tmp = sol.',field,';']);     
   else
      tmp = field;
   end

   % get the data depending on the viewing direction

   if (strcmp(opt.dir,'x'))
           slc = opt.xslc; c0 = sol.x;
      dir1 = 'y'; n1 = ny; c1 = sol.y;  % to be plotted on 1-x-axis of graph
      dir2 = 'z'; n2 = nz; c2 = sol.z;  % to be plotted on 2-y-axis of graph
   elseif (strcmp(opt.dir,'y'))
           slc = opt.yslc; c0 = sol.y;
      dir1 = 'x'; n1 = nx; c1 = sol.x;
      dir2 = 'z'; n2 = nz; c2 = sol.z;
   elseif (strcmp(opt.dir,'z'))
           slc = opt.zslc; c0 = sol.z;
      dir1 = 'x'; n1 = nx; c1 = sol.x;
      dir2 = 'y'; n2 = ny; c2 = sol.y;
   else
      disp(sprintf('ERROR: illegal value for opt.dir=%s',opt.dir));
      return;
   end

   % check size for requested coordinate direction
   if ((n1==1 | n2==1) & ...
       (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
      disp(sprintf('ERROR: the block has %d x %d points per %c-slice;', n1, n2, opt.dir));
      disp(sprintf('       %s-plot not available for %c-slices.', opt.typplot, opt.dir));
      return
   end

   % determine requested slice number
   if (isempty(slc))
      disp(sprintf('ERROR: no %c-coordinate provided for %c%c-slice',opt.dir,dir1,dir2));
      return;
   end
   if (ischar(slc))
      if (strcmp(lower(slc),'max'))
         islc = -1;
      else
         disp(sprintf('ERROR: unknown option slc=''%c''',slc));
         return;
      end
   else
      [mn, islc] = min(abs(c0-slc));
      if (abs(c0-slc)>1e-4);
         disp(sprintf('Showing results for %c=%5.3f, nearest to %cslc=%5.3f',...
                                                                        opt.dir, c0(islc), opt.dir, slc));
      end
   end

   % get the data for the slice
   if (islc<=0)
      idir = double(lower(opt.dir)) - double('x') + 1;
      dta = squeeze( max(abs(tmp), [], idir ));
   else
      if (strcmp(opt.dir,'x'))
         dta = reshape( tmp(islc,:,:), n1,n2);
      elseif (strcmp(opt.dir,'y'))
         dta = reshape( tmp(:,islc,:), n1,n2);
      elseif (strcmp(opt.dir,'z'))
         dta = reshape( tmp(:,:,islc), n1,n2);
      end
   end

   if (~isempty(opt.cntrlvl) && ~strcmp(opt.cntrlvl,'auto'))
      cntrlvl = opt.cntrlvl;
   end

   % prepare for logarithmic color scale
   if (strcmp(opt.scale,'log'))
      dta     = log10(abs(dta));
      if (~isempty(opt.cntrlvl) && ~strcmp(opt.cntrlvl,'auto'))
         ix      = find(cntrlvl>0);
         cntrlvl = log10(cntrlvl(ix));
      end
   end

   % for surf-plot: add one column and row to the data
   dta = [dta, dta(:,end)]; 
   dta = [dta; dta(end,:)];

   % avoid constant data
   if (max(max(dta))-min(min(dta))<1e-20)
      dta(1,1) = dta(1,1)+1e-20;
   end

   % interpret c1, c2 as cell-centers, compute corners of the cells for surf-plot
   if (length(c1)==1)
      d1 = (c2(end)-c2(1))/n2;
      cornr1 = [ c1(1)-d1(1)/2; c1(end)+d1(end)/2];
   else
      d1 = c1(2:end) - c1(1:end-1); % step size
      cornr1 = [ c1(1)-d1(1)/2; c1(1:end-1)+d1(1:end)/2; c1(end)+d1(end)/2];
   end
   if (length(c2)==1)
      d2 = (c1(end)-c1(1))/n1; 
      cornr2 = [ c2(1)-d2(1)/2; c2(end)+d2(end)/2];
   else
      d2 = c2(2:end) - c2(1:end-1); % step size
      cornr2 = [ c2(1)-d2(1)/2; c2(1:end-1)+d2(1:end)/2; c2(end)+d2(end)/2];
   end
   % do not extend with half a cell-size at z=0:
   if (strcmp(dir2,'z') & abs(c2(1))<1e-5), cornr2(1) = 0; end;
   if (strcmp(dir2,'z') & length(c2)>1 & abs(c2(end))<1e-5), cornr2(end) = 0; end;

   % plot contours or surf plot
   if (strcmp(opt.typplot,'contour'))
      if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
         [C, h] = contour(c1, c2, dta(1:end-1,1:end-1)');
      else
         [C, h] = contour(c1, c2, dta(1:end-1,1:end-1)', cntrlvl);
      end
   elseif (strcmp(opt.typplot,'contourf'))
      if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
         [C, h] = contourf(c1, c2, dta(1:end-1,1:end-1)');
      else
         [C, h] = contourf(c1, c2, dta(1:end-1,1:end-1)', cntrlvl);
      end
   else
      % size(cornr1), size(cornr2), size(dta)
      l1 = surf(cornr1, cornr2, dta');
      % shading faceted
      % set(l1,'linestyle','--');
      % l2=mesh(c1,c2,10*ones(size(dta)-1));
      % set(l2,'facecolor','none');
      % set(l2,'edgecolor','k');
      view([0 90]);
      axis([min(cornr1), max(cornr1), min(cornr2), max(cornr2)]);
   end

   if (isnumeric(opt.clabel) & (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
      clabel(C, h, opt.clabel);
   elseif (strcmp(opt.clabel,'on') & (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
      clabel(C, h);
   end

   if (~isempty(opt.colormap))
      colormap(opt.colormap);
   end

   xlabel(sprintf('%c_c [mm]',dir1));
   ylabel(sprintf('%c_c [mm]',dir2));
   % xlabel(sprintf('%c-coordinate [mm]',dir1));
   % ylabel(sprintf('%c-coordinate [mm]',dir2));
   if (islc>=1)
      title(sprintf(' at %c=%5.3f', opt.dir, c0(islc)));
   else
      title(sprintf(', max in %c', opt.dir));
   end

   % optionally adjust ticks

   % set(gca,'xtick', xcentr, 'xticklabel', num2str(xcentr', '%5.3f'));
   % set(gca,'ytick', ycentr, 'yticklabel', num2str(ycentr', '%5.3f'));

end % function show_2d_slice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = log_colorbar(hc, opt)
% 
% function [ ] = log_colorbar(hc, opt)
%
% change colorbar to logarithmic scale when necessary

   if (strcmp(opt.scale,'log'))
      if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
         yt=get(hc,'ytick');
      else
         ix = find(opt.cntrlvl>0);
         yt = log10(opt.cntrlvl(ix));
      end
      set(hc,'ytick',yt,'yticklabel',10.^yt);
   end
end % function log_colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % function plotstrs
