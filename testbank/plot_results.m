
if (~exist('loadcase'))
   addpath('../matlab');
end
if (~exist('plot_arrow'))
   addpath('../matlab_intern');
end

expnam = 'testbank';

% show_cases       = [ 18:23, 27:39 ]; % testbank.inp: cases using module 1
show_cases       = [ 35 ];
fig_ofs          = 0;
pause_after_plot = 0 * (length(show_cases)>1);
print_figures    = 0;

% list all cases used in testbank.inp

ncase  = 39;
module = 3 * ones(1, ncase);
nam    = repmat(' ', ncase,1);
ax1    = zeros(ncase, 4);
rfname = repmat(' ', ncase,1); rmirrory = zeros(1,ncase); rmirrorz = zeros(1,ncase);
wfname = repmat(' ', ncase,1); wmirrory = zeros(1,ncase); wmirrorz = zeros(1,ncase);

icase = 18; nam(icase,1:30)    = 'large roll & irregularities   '; module(icase)   = 1; 
            rfname(icase,1:30) = 'MBench_UIC60_v3.prr           '; rmirrory(icase) = 0;
            wfname(icase,1:30) = 'MBench_S1002_v3.prw           '; wmirrory(icase) = 0;
            ax1(icase,:) = [705 735 -15 5];
icase = 19; nam(icase,1:30)    = 'circular wheel/flat inclined  '; module(icase)   = 1;
            rfname(icase,1:30) = 'flat_plane.prr                '; rmirrory(icase) = 0;
            wfname(icase,1:30) = 'circ_r100.prw                 '; wmirrory(icase) = 0;
            ax1(icase,:) = [750 770 -15 5];
icase = 20; nam(icase,1:30)    = 'vertical section in wheel     '; module(icase)   = 1;
            rfname(icase,1:30) = 'os_l.csv                      '; rmirrory(icase) = 0;
            wfname(icase,1:30) = 'fr_ol.csv                     '; wmirrory(icase) = 0;
            ax1(icase,:) = [-35 -10 -12 6];
icase = 21; nam(icase,1:30)    = 'compressible sheet, combining '; module(icase)   = 1;
            rfname(icase,1:30) = 'meas_roller.prr               '; rmirrory(icase) = 0;
            wfname(icase,1:30) = 'final_clean.prw               '; wmirrory(icase) = 0;
icase = 22; nam(icase,1:30)    = 'steady curving example, high  '; module(icase)   = 1;
            rfname(icase,1:30) = 'avg_HR.BAN                    '; rmirrory(icase) = 1;
            wfname(icase,1:30) = 'avg_wheel.whl                 '; wmirrory(icase) = 0;
            ax1(icase,:) = [-730 -712 -3 12];
icase = 23; nam(icase,1:30)    = 'steady curving example, low   '; module(icase)   = 1;
            rfname(icase,1:30) = 'avg_low2.ban                  '; rmirrory(icase) = 1;
            wfname(icase,1:30) = 'avg_wheel.whl                 '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 740 760 -12 5];
icase = 27; nam(icase,1:30)    = 'iteration 288 with dents      '; module(icase)   = 1;
            rfname(icase,1:30) = 'rail_left_iter288.ban         '; rmirrory(icase) = 1;
            wfname(icase,1:30) = 'whl_left_iter288.whl          '; wmirrory(icase) = 0;
            ax1(icase,:) = [ -725 -705 8 26];
icase = 28; nam(icase,1:30)    = 'iteration 288, many kinks     '; module(icase)   = 1;
            rfname(icase,1:30) = 'rail_left_iter288.ban         '; rmirrory(icase) = 1;
            wfname(icase,1:30) = 'whl_left_iter288.whl          '; wmirrory(icase) = 0;
            ax1(icase,:) = [ -725 -705 8 26];
icase = 29; nam(icase,1:30)    = 'contact at bottom of flange   '; module(icase)   = 1;
            rfname(icase,1:30) = 'rail_left_iter795.ban         '; rmirrory(icase) = 1;
            wfname(icase,1:30) = 'whl_left_iter795.whl          '; wmirrory(icase) = 0;
            ax1(icase,:) = [ -715 -695 18 36];
icase = 30; nam(icase,1:30)    = 'guard rail, narrow pot.contact'; module(icase)   = 1;
            rfname(icase,1:30) = 'guard_rail_test.ban           '; rmirrory(icase) = 1; rmirrorz(icase) = 1;
            wfname(icase,1:30) = 'Car7813_0001r.whl             '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 650 780 -30 55];
icase = 31; nam(icase,1:30)    = 'guard, wheel with meas.noise/1'; module(icase)   = 1;
            rfname(icase,1:30) = 'guard_rail_test.ban           '; rmirrory(icase) = 1; rmirrorz(icase) = 1;
            wfname(icase,1:30) = 'Car7216_0001r.whl             '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 650 780 -30 55];
icase = 32; nam(icase,1:30)    = 'guard rail contact            '; module(icase)   = 1;
            rfname(icase,1:30) = 'guard_rail_test.ban           '; rmirrory(icase) = 1; rmirrorz(icase) = 1;
            wfname(icase,1:30) = 'Car7358_0001r.whl             '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 650 780 -30 55];
icase = 33; nam(icase,1:30)    = 'u.guide planar separation     '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 34; nam(icase,1:30)    = 'u.guide planar combination    '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 35; nam(icase,1:30)    = 'u.guide conformal combination '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 36; nam(icase,1:30)    = 'conformal with massless rail  '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 37; nam(icase,1:30)    = 'massless rail w. initial estim'; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 38; nam(icase,1:30)    = 'trans.shift, initial step     '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];
icase = 39; nam(icase,1:30)    = 'trans.shift, split patches    '; module(icase)   = 1;
            rfname(icase,1:30) = 'site_b_hr.prr                 '; rmirrory(icase) = 0; rmirrorz(icase) = 0;
            wfname(icase,1:30) = 'worn_high_mileage.prw         '; wmirrory(icase) = 0;
            ax1(icase,:) = [ 690 790 -20 40];

% determine version information on the Matlab graphics system

if (~verLessThan('matlab', '8.7.0'))
   set(0,'defaultlegendautoupdate', 'off');
end

for icase = show_cases

   if (module(icase)==1)
      % plot results for Module 1

      prr = read_profile(['profiles/', deblank(rfname(icase,:))], 0, rmirrory(icase), rmirrorz(icase));
      prw = read_profile(['profiles/', deblank(wfname(icase,:))], 1, wmirrory(icase), wmirrorz(icase));
      s = loadcase(expnam, icase, 0);

      opt = plot3d;
      opt.typplot  = 'rw_rear';
      opt.rw_surfc = 'both';
      opt.field    = 'pn';

      figure(fig_ofs+icase);
      plot3d(s, opt, prr, prw);
      title(sprintf('Case %d: %s', icase, deblank(nam(icase,:))));
      if (max(abs(ax1(icase,:)))>0), axis(ax1(icase,:)); end

   else
      % plot results for Module 3

      disp(sprintf('Case %d: module 3 n.y.a.',icase));

   end
   if (pause_after_plot), pause; end
end

clear ncase rfname rmirrory rmirrorz wfname wmirrory wmirrorz module ax1 nam;
clear pause_after_plot print_figures show_cases fig_ofs expnam icase s;

