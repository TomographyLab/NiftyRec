function f1 = readecatvol( fpname, matnums )
% vol = readecatvol( fpname, matnums )
%
% structure vol:
% vol.vol{ fr }( x, y, z )
% vol.hd{ fr }
% vol.sf( fr ) % does never contain decay correction
% vol.times( fr, 1:4 )
% vol.decaycf( fr ) % not always valid
% use matnums = [ -inf, -inf, -inf, -inf, -inf; inf, inf, inf, inf, inf ] 
%  in order to read the entire file.
% the columns in matnums refer to frame, plane, gate, data, bed


if size( matnums, 1 ) == 1
  matnums = [ matnums; matnums ];
end

f1 = [];
f1.fpname = fpname;
[ fid1, message ]       = ecatfile( 'open', f1.fpname );
[ mh1 ]                 = ecatfile( 'mh', fid1 );
[ matranges, message ] = ecatfile( 'matranges', fid1 );

if fid1 == -1
  f1 = [];
  return
end

matnums( 1, isinf( matnums( 1, : ) ) ) = ...
  matranges( 1, isinf( matnums( 1, : ) ) );

matnums( 2, isinf( matnums( 2, : ) ) ) = ...
  matranges( 2, isinf( matnums( 2, : ) ) );

f1.vol = {};
f1.hd = {};
for fr = matnums( 1, 1 ):matnums( 2, 1 ) %                     N O T E
  if strcmp( mh1.file_system, 'ecat6.4' )
    vol = [];
    for pl = matnums( 1, 2 ):matnums( 2, 2 )
      [ vol1, hd1, status ]    = ecatfile( 'read', fid1, [ fr pl 1 0 0 ] );
      if pl == matnums( 1, 2 )
        hd = hd1;
      end
      vol(:,:,pl) = vol1;
      sf1(pl) = hd1.sh.quant_scale;
    end
    % convert to common scale factor. Only one hd is kept.
    sf = max( sf1 );
    hd.sh.quant_scale = sf;
    hd.sh.plane_eff_corr_fctr = 1; %                         N O T E
    for pl = matnums( 1, 2 ):matnums( 2, 2 )
      vol(:,:,pl) = int16( round( double( round( vol( :, :, pl ) ) * ( sf1(pl) / sf ) ) ) );
    end
  else
    [ vol, hd, status ]    = ecatfile( 'read', fid1, [ fr 1 1 0 0 ] );
  end

   f1.vol{ fr - matnums( 1, 1 ) + 1 } = vol;
   f1.hd{ fr - matnums( 1, 1 ) + 1} = hd;
   %fprintf( '.' )
 end
status = ecatfile( 'close', fid1 );

f1.times = [];
% f1.sf = [];
f1.sf_nondecaycorr = [];
f1.decaycf = [];
for fr = 1:length(f1.hd)
    start = 0;%f1.hd{fr}.sh.frame_start_time;
    dur = 1;%f1.hd{fr}.sh.frame_duration;
    f1.times( end + 1, : ) = [ start, start + dur / 2, start + dur, dur ];
  if strcmp( mh1.file_system, 'ecat6.4' )
    % ecat6.4
    switch f1.hd{ fr }.mh.file_type
    case 2 % image
      f1.sf_nondecaycorr( end + 1 ) = f1.hd{fr}.sh.quant_scale / ...
        f1.hd{fr}.sh.decay_corr_fctr * ...
        f1.hd{fr}.sh.ecat_calibration_fctr;
      f1.decaycf( end + 1 ) = f1.hd{fr}.sh.decay_corr_fctr;
    end
  else
    % ecat7
    switch f1.hd{ fr }.mh.file_type
    case 11 % 3D sinogram
      f1.sf_nondecaycorr( end + 1 ) = f1.hd{fr}.sh.scale_factor;
      f1.decaycf( end + 1 ) = 1;
    case 7 % image
      f1.sf_nondecaycorr( end + 1 ) = f1.hd{fr}.sh.scale_factor / ...
        f1.hd{fr}.sh.decay_corr_fctr * ...
        f1.hd{fr}.mh.ecat_calibration_factor;
      f1.decaycf( end + 1 ) = f1.hd{fr}.sh.decay_corr_fctr;
    end
  end

end
f1.times = f1.times / 1000;
f1.sf = f1.sf_nondecaycorr;