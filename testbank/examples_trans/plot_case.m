prr = read_profile('../profiles/site_b_hr.prr');
prw = read_profile('../profiles/worn_high_mileage.prw');
s = loadcase('t3_t2_conformal',2);
opt = plot3d;
opt.typplot = 'rw_rear';
opt.rw_surfc = 'both';
opt.field = 'pn';
plot3d(s, opt, prr, prw);
