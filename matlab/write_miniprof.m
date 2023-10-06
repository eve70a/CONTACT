function [ ] = write_miniprof(p, fname, is_rail, use_yz, idebug)

% function [ ] = write_miniprof(p, fname, [is_rail], [use_yz], [idebug])
%
% write profile to file in Miniprof format 
%  - default (use_yz=0) : write data given in OrigData
%  - option  (use_yz>0) : overwrite OrigData with ProfileY and ProfileZ
%
% use e.g. p = struct('OrigData', [y, z], 'ColumnDef', 'X,Y');

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<3 | isempty(is_rail))
   is_rail = -1;
end
if (nargin<4 | isempty(use_yz))
   use_yz = 0;
end
if (nargin<5 | isempty(idebug))
   idebug = 0;
end

if (exist(fname,'file'))
   disp(['ERROR: file ',fname,' already exists, please remove']);
   return
end

% use ProfileY, ProfileZ if so requested, else OrigData will be used

if (use_yz)
   p.OrigData = [ p.ProfileY, p.ProfileZ ];
   p.ColumnDef = 'X,Y';
end

f = fopen(fname, 'w');

known_fields = struct( ...
   'ProfileData',                   'internal', ...
   'ProfileY',                      'internal', ...
   'ProfileZ',                      'internal', ...
   'ProfileS',                      'internal', ...
   'ProfileAngle',                  'internal', ...
   'ProfileCurvature',              'internal', ...
   'ProfileKYield',                 'internal', ...
   'OrigData',                      'internal', ...
   'Mask',                          'internal', ...
   'Fname',                         'internal', ...
   'FileName',                      'str', ...
   'ProgramName',                   'str', ...
   'ProgramDate',                   'str', ...
   'ProgramVer',                    'str', ...
   'Date',                          'str', ...
   'Time',                          'str', ...
   'MPCalDate',                     'str', ...
   'MPCalTime',                     'str', ...
   'UserName',                      'str', ...
   'EmplNum',                       'int', ...
   'Track',                         'str', ...
   'Line',                          'str', ...
   'Direction',                     'str', ...
   'Chainage',                      'str', ...
   'Curve',                         'str', ...
   'Position',                      'str', ...
   'Location_Latitude',             'real', ...
   'Location_Longitude',            'real', ...
   'Location_Altitude',             'real', ...
   'Location_HorizontalAccuracy',   'real', ...
   'Location_VerticalAccuracy',     'real', ...
   'Rail',                          'str', ...
   'TrackSection',                  'unknown', ...
   'KP',                            'unknown', ...
   'Stock',                         'str', ...
   'CarNo',                         'int', ...
   'Bogie',                         'int', ...
   'SurfaceCondition',              'unknown', ...
   'Side',                          'str', ...
   'AxleNo',                        'int', ...
   'Mileage',                       'str', ...
   'WheelID',                       'int', ...
   'RCF',                           'unknown', ...
   'Flat',                          'unknown', ...
   'ReferenceProfile',              'str', ...
   'ProfileAlignment',              'unknown', ...
   'Instrument_Firmware',           'str', ...
   'Instrument_Connection',         'int', ...
   'Windows_Version',               'str', ...
   'Windows_ServicePack',           'str', ...
   'Windows_NetFramework',          'str', ...
   'Windows_ComputerName',          'str', ...
   'Windows_UserName',              'str', ...
   'Battery_Level',                 'real', ...
   'Battery_Capacity',              'real', ...
   'Comment',                       'str', ...
   'Flag',                          'unknown', ...
   'MPTypeNo',                      'int', ...
   'MPSerNo',                       'int', ...
   'Xoffset',                       'real', ...
   'Yoffset',                       'real', ...
   'Elevation',                     'real', ...
   'Tilt',                          'unknown', ...
   'AreaRef',                       'unknown', ...
   'AreaGain',                      'unknown', ...
   'AreaLoss',                      'unknown', ...
   'BaseTilt',                      'unknown', ...
   'BaseOffset',                    'unknown', ...
   'RefPoint1',                     'real', ...
   'RefPoint2',                     'real', ...
   'RefPoint3',                     'real', ...
   'RefPoint4',                     'real', ...
   'RefPoint5',                     'real', ...
   'RefPoint6',                     'real', ...
   'RefPoint7',                     'real', ...
   'RefPoint8',                     'real', ...
   'RefPoint9',                     'real', ...
   'SectionNo',                     'unknown', ...
   'CrownR',                        'real', ...
   'CrownRadius',                   'real', ...
   'TransformState',                'str', ...
   'Transformation',                'str', ...
   'ProfileShaping',                'unknown', ...
   'OriginalHeader',                'unknown', ...
   'ColumnDef',                     'str', ...
   'XYPoints',                      'int', ...
   'RailAlignType',                 'int', ...
   'RailAngle',                     'real', ...
   'RailAngle_AlarmStatus',         'int', ...
   'AlignDX',                       'unknown', ...
   'AlignDY',                       'unknown', ...
   'AlignDA',                       'unknown', ...
   'AlignOX',                       'unknown', ...
   'AlignOY',                       'unknown', ...
   'W1',                            'real', ...
   'W2',                            'real', ...
   'W3',                            'real', ...
   'W1_AlarmLevel',                 'real', ...
   'W1_AlarmStatus',                'int', ...
   'W2_AlarmStatus',                'int', ...
   'W3_AlarmStatus',                'int', ...
   'UOF1',                          'str', ...
   'UOF2',                          'str', ...
   'UOF3',                          'str', ...
   'UOF4',                          'str', ...
   'Chars',                         'str', ...
   'Temperature',                   'real', ...
   'GradeRatio',                    'real', ...
   'Grade',                         'real', ...
   'SuperElevationHeight',          'real', ...
   'SuperElevation',                'real', ...
   'RailGaugePoint',                'str', ...
   'RailRodLength',                 'str', ...
   'RodLength',                     'real', ...
   'GaugeRodLength',                'real', ...
   'GaugeSensorDistance',           'real', ...
   'GaugeSensorStdDev',             'real', ...
   'GaugeSensorSamples',            'real', ...
   'WheelAdjustPoint',              'real', ...
   'WheelDiameterTaperline',        'real', ...
   'AlarmWarningCount',             'int', ...
   'AlarmFailureCount',             'int', ...
   'Gauge',                         'real', ...
   'Gauge_AlarmStatus',             'int', ...
   'DiameterFlange',                'real', ...
   'DiameterTaperline',             'real', ...
   'DiameterFlange_AlarmStatus',    'int', ...
   'DiameterTaperline_AlarmStatus', 'int', ...
   'Sd',                            'real', ...
   'Sh',                            'real', ... 
   'qR',                            'real', ...
   'Sd_AlarmStatus',                'int', ...
   'Sh_AlarmStatus',                'int', ... 
   'qR_AlarmStatus',                'int'  ... 
);

% first write all the metadata to the file: keyword=value

fields_p = fieldnames(p);
num_fields = length(fields_p);

for ifld = 1 : num_fields

   my_field_name = check_field(known_fields, fields_p(ifld));

   if (~isempty(my_field_name))
      % interpret numerical value
      if (strcmp(known_fields.(my_field_name), 'internal'))
         % disp(sprintf('Generated field "%s"', fields_p{ifld}));
      elseif (strcmp(known_fields.(my_field_name), 'int'))
         nval = length(p.(my_field_name));
         fmt = '%s=%d';
         for ival = 2 : nval
            fmt = [fmt, ', %d'];
         end
         fmt = [fmt, '\n'];
         fprintf(f, fmt, my_field_name, p.(my_field_name) );
      elseif (strcmp(known_fields.(my_field_name), 'real') | ...
              strcmp(known_fields.(my_field_name), 'num'))
         % TODO: arrays!
         nval = length(p.(my_field_name));
         fmt = '%s=%6.4f';
         for ival = 2 : nval
            fmt = [fmt, ', %6.4f'];
         end
         fmt = [fmt, '\n'];
         fprintf(f, fmt, my_field_name, p.(my_field_name) );
      % store string value
      elseif (strcmp(known_fields.(my_field_name), 'str') | ...
              strcmp(known_fields.(my_field_name), 'unknown'))
         fprintf(f, '%s=%s\n', my_field_name, p.(my_field_name) );
      % unknown type: error
      else
         disp(sprintf('Field %s has unknown field-type %s.',my_field_name, ...
                known_fields.(my_field_name) ));
      end
   % unknown keyword: report
   else
      disp(sprintf('Unknown field "%s"',fields_p{ifld}));
   end
end % all fields

% create format on basis of ColumnDef

coldef = upper(p.ColumnDef);

fmt_X = ' %8.4f';
fmt_Y = ' %8.4f';
fmt_A = ' %7.4f';
fmt_C = ' %18.15f';
fmt_N = ' %3d';
fmt_K = ' %7.3f';
str   = ['fmt_', strrep(coldef, ',', ', fmt_') ];
fmt   = eval([ '[',str,']' ]);
fmt   = [fmt, ' \n'];

% Write the profiledata to the file

fprintf(f, '\n');
for ipoint = 1 : size( p.OrigData,1 )
   fprintf(f, fmt, p.OrigData(ipoint,:));
end

fclose(f);


end % function write_miniprof

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ my_field ] = check_field( known_fields, field_in )

   fn = fieldnames(known_fields);
   chk = strcmpi(field_in, fn);
   i = find(chk);
   if (isempty(i))
      my_field = [];
   else
      my_field = fn{i};
   end

end % function check_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

