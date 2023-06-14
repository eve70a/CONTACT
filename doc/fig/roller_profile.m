
function [xprf, yprf] = roller_profile( r_nom, r_axle );

% function [xprf, yprf] = roller_profile( r_nom, r_axle );

% points obtained from ginput:

x = [ 253.6670 221.7798 200.9555 192.4957 185.3373 182.0835 183.3850 187.9403 225.6844 230.8905 ...
      233.4935 234.1443 230.8905 219.8275 185.0000 180.0000 180.0000 325.0000 325.0000 320.0000 ...
      286.8557 279.6974 275.1421 273.1898 277.0944 280.9989 321.3460 323.9490 325.2505 325.2505 ...
      320.0445 312.2354 297.2679 253.6670 ];
y = [ 337.7863 338.4371 333.2310 328.6757 315.6605 272.7104 267.5043 262.9490 250.5846 245.3785 ...
      238.8709  93.1009  80.0857  69.0228  50.0000  42.0000  20.0000  20.0000  42.0000  50.0000 ...
       70.3243  78.1334  89.1963 239.5217 245.3785 249.2831 264.2505 265.5521 270.1074 311.1052 ...
      320.8666 329.3265 333.8818 337.7863 ];

if (0==1)
   figure(2); clf; hold on;
   plot( x,  y, '-*')
   % plot(-x,  y, '-o')
   grid on
   for i=1:length(x),
      text(x(i), y(i), num2str(i));
   end
end

% shift center to x=0

x = x - 254;

% create symmetric profile

xprf = (x - fliplr(x))/2;
yprf = (y + fliplr(y))/2;

% set y positive downwards, origin at top of rail

yprf = - (yprf - max(yprf));

% scale such that foot is x= -70:70 at y=160, 

xprf = xprf *  70/max(xprf);
yprf = yprf * 160/max(yprf);

% extend the web, such that the foot is at y=r_nom - r_axle

ix_foot = find(yprf>100);
yprf(ix_foot) = yprf(ix_foot) - 160 + r_nom - r_axle;

if (0==1)
   figure(3); clf; hold on;
   % plot(-x,  y, '-o')
   plot(xprf, yprf, '-o')
   grid on
end
