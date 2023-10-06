
print_fig = 1;

fkin = 0.15;
fslp = 0.15;
shlf = 3;
s  = [0:0.1:14];
ix = [1:10:141];

f1 = fkin + fslp * ones(size(s));
f2 = fkin + fslp * max(0, 1-s/(2*shlf));
f3 = fkin + fslp ./ (1 + s/shlf);
f4 = fkin + fslp ./ (1 + (s/shlf).^2);
f5 = fkin + fslp * exp(-log(2)*s/shlf);

figure(1); clf; hold on;

l = plot(s, [f1; f2; f4; f5; f3]);
set(l(3), 'linestyle', '--');
set(l(5), 'linestyle', '-.');

xlabel('Absolute slip velocity s_a [m/s]');
ylabel('Coefficient of friction \mu(s_a) [-]');
axis([0 14 0.13 0.32]);
grid on;

set(gca,'colororderindex',1);
plot(s(ix), f1(ix), '+');
plot(s(ix), f2(ix), 'd');
plot(s(ix), f4(ix), 'x');
plot(s(ix), f5(ix), '*');
plot(s(ix), f3(ix), 'o');

[lg,ln] = legend(l([1,2,5,3,4]), 'Coulomb friction', 'Linear+constant', ...
        'Rational law 1', 'Rational law 2', 'Exponential law', ...
        'location','northeast');

set(ln(7), 'marker','+');
set(ln(9), 'marker','d');
set(ln(11), 'marker','o');
set(ln(13), 'marker','x');
set(ln(15), 'marker','*');

text(0.95, 0.625, '\mu_{kin} = 0.15, \mu_1 = 0.15, s_{hlf} = 3 m/s', ...
                'units','normalized', 'horizontalalignment','right');

if (print_fig)
   print -djpeg95 fric_laws.jpg
end

