function [ res1, res2, res3 ] = ecatfile( action, p1, p2, p3, p4, p5 )
% USAGE:
% [ fid, message ]       = ecatfile( 'open', filename, file_system );
% [ fid, message ]       = ecatfile( 'open', filename, 'ecat7' );
% [ fid, message ]       = ecatfile( 'open', filename, 'ecat6.4' );
% [ fid, message ]       = ecatfile( 'open', filename, '' );
% [ fid, message ]       = ecatfile( 'create', filename, mh );
% [ mh, message ]        = ecatfile( 'mh', fid );
% [ vol, hd, message ]   = ecatfile( 'read', fid, selmatrix );
% [ vol, hd, message ]   = ecatfile( 'read', fid, selmatrix, segment ); % 3D sinogram
% [ descr, hd, message ] = ecatfile( 'read_offset', fid, selmatrix );
% [ descr, hd, message ] = ecatfile( 'read_offset', fid, selmatrix, segment ); % 3D sinogram
% [ hd, message ]        = ecatfile( 'hd', fid, selmatrix );
% message                = ecatfile( 'writemh', fid, mh );
% message                = ecatfile( 'write', fid, vol, hd, selmatrix );
% message                = ecatfile( 'write', fid, vol, hd, selmatrix, segment ); % 3D sinogram
% message                = ecatfile( 'close', fid );
% message                = ecatfile( 'close', 'all' );
% [ matranges, message ] = ecatfile( 'matranges', fid );
% [ matlist, matstatus, message ] = ecatfile( 'matlist', fid );
% [ hd, hds, message ]   = ecatfile( 'blankhd', file_system );
% [ mh, hds, message ]   = ecatfile( 'blankhd', 'ecat7' );
% [ mh, hds, message ]   = ecatfile( 'blankhd', 'ecat6.4' );
% [ sh, hds, message ]   = ecatfile( 'blankhd', 'ecat7', file_type );
% [ sh, hds, message ]   = ecatfile( 'blankhd', 'ecat6.4', file_type );
%
% The user has the responsibility that the written headers are compatible
%   with the file and data structure
% Revision: April 24, 2001.
% ecatfile has been tested under MATLAB5.3. It has not been tested thoroughly under MATLAB6.
% 1999, 2000, 2001 by Flemming Hermansen.

% Known problem with opened files:
% The ecatfile.m reads all headers and all the index blocks when the ECAT file is opened. 
% Further calls to open will not reread the headers and index blocks if the file is 
% already open. Each call to open or create will generate different fid values for the 
% same file. The headers and the index blocks will only be discarded, when all fid's have 
% been closed (or when 'close all' is issued). This will give problems if the file is 
% modified outside the Matlab process, because the headers and index blocks in the Maltab 
% process then don't correspond to those in the physical file. This may also happen, if a 
% file is renamed or copied. In this situation the user may call 
%    message = ecatfile( 'close', 'all' )
% to make sure that ecatfile.m will reread the headers and index blocks. Usually, each 
% program should close a file after using it. Still, files may be left opened, if a 
% program crashes before closing the file. A subroutine that is used to read portions of a 
% file will thus open and close the file many times. This will impose a high overhead 
% because the headers and index blocks will be read repeatedly. This overhead may be 
% prevented if the calling program, itself, opens the file before calling the 
% subroutine, and closes the file when it has performed the last call to the subroutine.

% UNRESOLVED PROBLEM: I have not seen a definite report from CTI stating all possible
% formats. There may, therefore, be errors due to misunderstandings and ignorance of 
% certain formats. 
% Especially, I do not know the machine format of the file. I have assumed that 6.4
% files can be read by using machineformat 'vaxg' in the fopen routine. Likewise, I
% assume that ecat7 files can be opened using 'ieee-be'. I assume that the entire file
% be read by the same machineformat.

% The routines have been tested, and they can read and write the
% following ecat6.4 files types: 1:*.scn, 2:*.img, 3:*.atn, 4:*.nrm
% and the following ecat7 files types: 3:*.a, 7:*.v, 11:*.S, 14:*.S
% It can read the header from the following file type: 13:*.N
% It should also be able to read other headers (see the tables below), but they have not been tested.
%
% The written files have been tested binarily, and there are some differences:
% The CTI generated files do not contain a truly doubly linked list of index blocks, as index block 2
% points back to block 0.
% Floating point numbers are not reproduced binarily exact, but the read value is the same, 
% except for values very close to zero. (assumedly a problem with the representation of 
% underflow, or illegal values)
% The last block in CTI generated ecat7 3D sinograms are filled with non-zeros.
% The last block in CTI generated ecat7 3D attenuation is truncated.

%NOT IMPLEMENTED and known errors:
% Problem not yet solved (requires a call to the c-routine stat in stat.h):
% A file may be opened twice or more times at the same time. The filename must, 
% however, be spelled exactly the same way. Otherwise, the ecatfile-system 
% can not keep the headers and directoryblocks updated for all fids.
% 
% check that seek succeeds.
% seek during write should write zeros in front of the matrix, if the file is too short.
%
% The system updates the index blocks after the data have been written, so that
% the file always will contain valid data. The only exception is when 3D files are written.
% Then the index block is updated after the first 3D-segment has been written. The file
% will thus be too short, until the last 3D-segment has been written.

% ------------------------------------------------------

% pers.handlelist
% pers.filelist

% pers.filelist{}.filename
% pers.filelist{}.mainheader
% pers.filelist{}.dirblocks.data( 1:n )
% pers.filelist{}.dirblocks.filepos( 1:n )
% pers.filelist{}.pos
% pers.filelist{}.pos.dirblock()
% pers.filelist{}.pos.fpgdb()
% pers.filelist{}.pos.blockstart()
% pers.filelist{}.pos.blockend()
% pers.filelist{}.pos.status()
% pers.filelist{}.pos.header{} % may be empty

   global pers
%   persistent pers

if isempty( pers )
   pers = [];
   pers.handlelist = [];
   pers.filelist = [];
end

switch action
   
case 'blankhd'
% [ hd, hds, message ]   = ecatfile( 'blankhd', file_system );
% [ mh, hds, message ]   = ecatfile( 'blankhd', 'ecat7' );
% [ mh, hds, message ]   = ecatfile( 'blankhd', 'ecat6.4' );
   file_system = p1;
   if nargin == 2
      [ hds , message ] = getmhs( file_system );
      if isempty( message )
         [ hd , message ] = blankhd( hds );
      end
      res1 = hd;
      res2 = hds;
      res3 = message;
      return
   end
   file_type = p2;
   [ hds , message ] = getshs( file_system, file_type );
   if isempty( message )
      [ hd , message ] = blankhd( hds );
   end
   res1 = hd;
   res2 = hds;
   res3 = message;
   
case 'matlist'
   % [ matlist, matstatus, message ] = ecatfile( 'matlist', fid );
   res1 = [];
   res2 = [];
   fid = p1;
   
   res3 = 'Illegal handle';
   if fid > length( pers.handlelist ) | fid < 1
      res3 = 'Illegal handle';
      return;
   elseif pers.handlelist( fid ) == 0
      res3 = 'Illegal handle';
      return;
   end
   file = pers.filelist{ pers.handlelist( fid )  };
   
   res1 = file.pos.fpgdb;
   res2 = file.pos.status;
   res3 = '';
   return
   
case 'matranges'
   % [ matranges, message ] = ecatfile( 'matranges', fid );
   res1 = [];
   res2 = '';
   fid = p1;
   
   if fid > length( pers.handlelist ) | fid < 1
      res2 = 'Illegal handle';
      return;
   elseif pers.handlelist( fid ) == 0
      res2 = 'Illegal handle';
      return;
   end
   file = pers.filelist{ pers.handlelist( fid )  };
   
   [ ranges, message ] = getranges( file );
   res1 = ranges;
   res2 = message;
   return
   
case 'close'
   % message                 = ecatfile( 'close', fid );
   % message                 = ecatfile( 'close', 'all' );
   
   res1 = '';
   fid = p1;
   if strcmp( fid,'all' )
      pers = [];
      pers.handlelist = [];
      pers.filelist = [];
      res1 = '';
      return
   end
   
   if fid > length( pers.handlelist ) | fid < 1
      res1 = 'Illegal handle';
      return;
   elseif pers.handlelist( fid ) == 0
      res1 = 'Illegal handle';
      return;
   end
   
   fileno = pers.handlelist( fid );
   pers.handlelist( fid ) = 0;
   
   if any( pers.handlelist == fileno ) | fid < 1
      return
   end
   
   pers.filelist{ fileno } = [];
   return
   
case 'write'
   % message = ecatfile( 'write', fid, vol, hd, selmatrix );
   % message = ecatfile( 'write', fid, vol, hd, selmatrix, segment ); % 3D sinogram
   fid = p1;
   vol = p2;
   hd = p3;
   selmatrix = p4;
   % 3D sinogram segment:
   if nargin < 6
      selsegment = 1;
   else
      selsegment = p5;
   end

   res1 = '';
   
   postponeindexupdate = 0;
   
   xyz_dimension = size( vol );
   while length( xyz_dimension ) < 3
      xyz_dimension = [ xyz_dimension, 1 ];
   end
   if prod( xyz_dimension ) ~= prod( xyz_dimension(1:3 ) )
      error( 'dimension of vol is higher than 3' )
      return
   end
   
   if fid > length( pers.handlelist )
      res1 = 'Illegal handle';
      return;
   elseif pers.handlelist( fid ) == 0
      res1 = 'Illegal handle';
      return;
   end
   
   file = pers.filelist{ pers.handlelist( fid ) };
   
   [ datadescr, message ] = getdatadescr( file.mainheader.file_system, ...
      file.mainheader.file_type, hd.sh.data_type );
   
   if isempty( message )
      [ hd, offsetfromheader, writefill, numdatablocks, message ] = ...
         set_xyz_dimension( hd, datadescr, xyz_dimension, selsegment  );
   end
   
   if isempty( message )
      [ hds, message ] = getshs( file.mainheader.file_system, file.mainheader.file_type );
   end
   
   if ~isempty( message )
      return
   end
   
   % total number of blocks == numdatablocks + hds.header_size / 512
   
   fidphys = fopen( file.filename, 'r+', file.mhmachineformat );
   
   if isempty( file.pos.fpgdb )
      matpos = 0;
   else
      matpos = firstno( ...
         file.pos.fpgdb( 1, : ) == selmatrix( 1 ) & ...
         file.pos.fpgdb( 2, : ) == selmatrix( 2 ) & ...
         file.pos.fpgdb( 3, : ) == selmatrix( 3 ) & ...
         file.pos.fpgdb( 4, : ) == selmatrix( 4 ) & ...
         file.pos.fpgdb( 5, : ) == selmatrix( 5 ) );
   end
   
   % calculate number of blocks
   %   data_type.ftype = 'int16';
   %   data_type.size = 2;
   %   data_type.fopentype = 'ieee-be';
   
   if matpos == 0
      % entry does not exist, create entry
      % allocate space in the index block, 
      % do not write the updated index block yet
      % function matnum = fpgdb2matnum( fpgdb )
      
      % find an index block that is not full
      indexblockno = firstno( file.dirblocks.data( 1, : ) ~= 0 );
      
      if indexblockno == 0
         % all index blocks are full, create a new index block
         newindexblock = max( max( file.dirblocks.data( 1, : ) ), ...
            max( file.pos.blockend ) ) + 1;
         lastindexblock = file.dirblocks.filepos( end );
         file.dirblocks.data( 2, end ) = newindexblock;
         file.dirblocks.data( 3, 1 ) = newindexblock;
         file.dirblocks.filepos = [ file.dirblocks.filepos, newindexblock ];
         file.dirblocks.data( :, end + 1 ) = ...
            [ 31, 2, lastindexblock, 0, zeros( 1, 124 ) ];
         
         % the following writing sequence ensures, that next will be consistent in 
         % the links, even in case of interruption during writing
         % write new index block;
         status = fseek( fidphys, ( newindexblock - 1 ) * 512, -1 );
         if status ~= 0
            fclose( fidphys );
            res1 = 'Error writing the index block';
            return
         end
         count = fwrite( fidphys, file.dirblocks.data( :, end ), 'uint32' );
         if count ~= 128
            fclose( fidphys );
            res1 = 'Error writing the index block';
            return
         end
         % write first index block
         status = fseek( fidphys, ( 2 - 1 ) * 512, -1 );
         if status ~= 0
            fclose( fidphys );
            res1 = 'Error writing the index block';
            return
         end
         count = fwrite( fidphys, file.dirblocks.data( :, 1 ), 'uint32' );
         if count ~= 128
            res1 = 'Error writing the index block';
            return
         end
         % write last index block
         if length( file.dirblocks.filepos ) > 2
            status = fseek( fidphys, ( file.dirblocks.filepos( end - 1 ) - 1 ) * 512, -1 );
            if status ~= 0
               fclose( fidphys );
               res1 = 'Error writing the index block';
               return
            end
            count = fwrite( fidphys, file.dirblocks.data( :, end - 1 ), ...
               'uint32' );
            if count ~= 128
               fclose( fidphys );
               res1 = 'Error writing the index block';
               return
            end
         end
         pers.filelist{ pers.handlelist( fid )  } = file;
         
         indexblockno = length( file.dirblocks.filepos );
      end % create a new index block
      
      % add item into indexblock
      % 0: matnum
      % 1: first data block (subheader)
      % 2: last data block
      % 3: status: 1:exists; 0: allocated, no data yet written, -1: matrix deleted
      startblock = max( max( file.dirblocks.data( 2, : ) ), ...
         max( [ file.pos.blockend, 0 ] ) ) + 1;
      p = file.dirblocks.data( 4, indexblockno ) * 4 + 4;
      file.dirblocks.data( p+1:p+4, indexblockno ) = [ ...
            fpgdb2matnum( selmatrix ), ...
            startblock, ...
            startblock + numdatablocks + hds.header_size / 512 - 1, ...
            1 ];
      file.dirblocks.data( 1, indexblockno ) = file.dirblocks.data( 1, indexblockno ) - 1;
      file.dirblocks.data( 4, indexblockno ) = file.dirblocks.data( 4, indexblockno ) + 1;
      
      matpos = length( file.pos.dirblock ) + 1;
      file.pos.dirblock( matpos ) = indexblockno;
      file.pos.fpgdb( :, matpos ) = selmatrix( : );
      file.pos.blockstart( matpos ) = startblock;
      file.pos.blockend( matpos ) = ...
         startblock + numdatablocks + hds.header_size / 512 - 1;
      file.pos.status( matpos ) = 1;
      file.pos.header{ matpos } = hd.sh;
      file.pos.relpos( matpos ) = p / 4;
      
      postponeindexupdate = 1;
      % postpone the updating of file and writing the above changes to the hard disk
      % until the matrix has been succesfully written
   else
      % entry exists, check size of existing entry
      if file.pos.blockend( matpos ) - file.pos.blockstart( matpos ) + 1 ~= ...
            numdatablocks + hds.header_size / 512
         error( 'Trying to change the size of a matrix' )
         return
      end      
      file.pos.header{ matpos } = hd.sh;
      if file.pos.status ~= 1
         % change status to 1
         p = file.pos.relpos( matpos );
         file.pos.status( matpos ) = 1;
         file.dirblocks.data( p * 4 + 4, indexblockno ) = 1;
         postponeindexupdate = 1;
      end
   end
   
   % pers.filelist{}.pos.dirblock()
   % pers.filelist{}.pos.fpgdb()
   % pers.filelist{}.pos.blockstart()
   % pers.filelist{}.pos.blockend()
   % pers.filelist{}.pos.status()
   % pers.filelist{}.pos.header{} % may be empty
   
   status = fseek( fidphys, ( file.pos.blockstart( matpos ) - 1 ) * 512, -1 );
   if status ~= 0
      fclose( fidphys );
      res1 = 'Error writing matrix';
      return
   end
   
   message = writehd( fidphys, hd.sh, hds );
   if ~isempty( message )
      fclose( fidphys );
      res1 = message;
      return
   end
   
   % 'a', ftell( fidphys )
   
   % skip to the right 3D segment
   if offsetfromheader ~= 0
      fseek( fidphys, offsetfromheader, 0 );
   end
   
   for pl = 1:xyz_dimension( 3 )
      count = fwrite( fidphys, double( vol( :, :, pl ) ), datadescr.ftype );
      if count ~= prod( xyz_dimension( 1:2 ) )
         fclose( fidphys );
         res1 = 'Error writing matrix';
         return
      end
   end
   
   % fill last block with zeros
   count = fwrite( fidphys, zeros( 1, writefill / datadescr.size ), datadescr.ftype );
   if count ~= writefill / datadescr.size
      fclose( fidphys );
      res1 = 'Error writing matrix';
      return
   end
   
   if hds.header_size + ...
         ( file.pos.blockstart( matpos ) - 1 + numdatablocks ) * 512 ...
         ~= ftell( fidphys )
      % warning( 'The length of the written matrix is wrong.' );
   end
   
   % write the updated indexblock
   if postponeindexupdate
      status = fseek( fidphys, ( file.dirblocks.filepos( indexblockno ) - 1 ) * 512, -1 );
      if status ~= 0
         fclose( fidphys );
         res1 = 'Error writing the index block';
         return
      end
      count = fwrite( fidphys, file.dirblocks.data( :, indexblockno ), 'uint32' );
      if count ~= 128
         fclose( fidphys );
         res1 = 'Error writing the index block';
         return
      end
   end
   
   pers.filelist{ pers.handlelist( fid )  } = file;
   
   fclose( fidphys );
   return
case { 'read', 'hd', 'read_offset' }
% [ vol, hd, message ]   = ecatfile( 'read', fid, selmatrix );
% [ vol, hd, message ]   = ecatfile( 'read', fid, selmatrix, segment ); % 3D sinogram
% [ descr, hd, message ] = ecatfile( 'read_offset', fid, selmatrix );
% [ descr, hd, message ] = ecatfile( 'read_offset', fid, selmatrix, segment ); % 3D sinogram
% [ hd, message ]        = ecatfile( 'hd', fid, selmatrix );
   fid = p1;
   selmatrix = p2;
   % 3D sinogram segment:
   if nargin < 4
      selsegment = 1;
   else
      selsegment = p3;
   end
   vol = [ ];
   hd = [ ];
   descr = [];
   message = '';
   if fid > length( pers.handlelist ) | fid <= 0
      message = 'Illegal handle';
   elseif pers.handlelist( fid ) == 0
      message = 'Illegal handle';
   else
      file = pers.filelist{ pers.handlelist( fid )  };
      if isempty( file.pos.fpgdb )
         message = 'Matrix does not exist';
      else
         matpos = firstno( ...
            file.pos.fpgdb( 1, : ) == selmatrix( 1 ) & ...
            file.pos.fpgdb( 2, : ) == selmatrix( 2 ) & ...
            file.pos.fpgdb( 3, : ) == selmatrix( 3 ) & ...
            file.pos.fpgdb( 4, : ) == selmatrix( 4 ) & ...
            file.pos.fpgdb( 5, : ) == selmatrix( 5 ) );
         if matpos == 0
            message = 'Matrix does not exist';
         end
      end
   end
   
   % pers.filelist{}.pos.dirblock()
   % pers.filelist{}.pos.fpgdb()
   % pers.filelist{}.pos.blockstart()
   % pers.filelist{}.pos.blockend()
   % pers.filelist{}.pos.status()
   % pers.filelist{}.pos.header{} % may be empty
   
   if isempty( message )
      fidphys = fopen( file.filename, 'r', file.mhmachineformat );
      fseek( fidphys, ( file.pos.blockstart( matpos ) - 1 ) * 512, -1 );
      [ hds, message ] = getshs( file.mainheader.file_system, file.mainheader.file_type );
   else
      fidphys = -1;
   end
   
   if isempty( message )
      [ sh, message ] = readhd( fidphys, hds );
      hd.sh =sh;
      hd.mh = file.mainheader;
   end
   
   if strcmp( action, 'hd' )
      res1 = hd;
      res2 = message;
      if fidphys > 0
         fclose( fidphys );
      end
      return
   end
   
   % 'a', ftell( fidphys )
   
   if isempty( message )
      [ descr, message ] = getdatadescr( file.mainheader.file_system, file.mainheader.file_type, sh.data_type );
   end
   
   if isempty( message )
      [ xyz_dimension, offsetfromheader ] = get_xyz_dimension( hd, descr, selsegment );
      xyz_dimension = [288   144    63];
      descr.xyz_dimension = xyz_dimension;
      descr.headeroffset = ( file.pos.blockstart( matpos ) - 1 ) * 512;
      descr.header_size = hds.header_size;
      descr.dataoffset = descr.headeroffset + hds.header_size + offsetfromheader;
   end
   
   if strcmp( action, 'read' ) & isempty( message )
      
      % skip to the right 3D segment
      if offsetfromheader ~= 0
         fseek( fidphys, offsetfromheader, 0 );
      end
      
      % preallocate vol:
      switch descr.ftype
      case 'float32'
         vol1 = zeros( xyz_dimension( 1 ), xyz_dimension( 2 ) );
      case { 'int8', 'int16' }
         vol1 = int16( zeros( xyz_dimension( 1 ), xyz_dimension( 2 ) ) );
      end
      vol = vol1( :, :, ones( xyz_dimension( 3 ), 1 ) );
      for pl = 1:xyz_dimension( 3 )
         [ vol1, count ] = fread( fidphys, prod( xyz_dimension( 1:2 ) ), ...
            descr.ftype );
         if count ~= prod( xyz_dimension( 1:2 ) ) 
            message = 'file too short';
            break
         end
         switch descr.ftype
         case 'float32'
            vol( :, :, pl ) = ...
               reshape( vol1, xyz_dimension( 1:2 ) );
         case { 'int8', 'int16' }
            vol( :, :, pl ) = ...
               reshape( int16( vol1 ), xyz_dimension( 1:2 ) );
         end
      end
   end 
   
   if strcmp( action, 'read' )
      if ~isempty( message )
         vol = [];
      end
      res1 = vol;
   else
      res1 = descr;
   end
   
   res2 = hd;
   res3 = message;
   if fidphys > 0
      fclose( fidphys );
   end
   
   return
   
case 'writemh'
   % message = ecatfile( 'writemh', fidphys, mh );
   fid = p1;
   mh = p2;
   res1 = '';
   if fid > length( pers.handlelist ) | fid <= 0
      res1 = 'Illegal handle';
   elseif pers.handlelist( fid ) == 0
      res1 = 'Illegal handle';
   end
   file = pers.filelist{ pers.handlelist( fid ) };
   fidphys = fopen( file.filename, ...
      'r+', file.mhmachineformat );
   if fidphys <= 0
      res1 = [ 'Can''t open file: ', file.filename ];
      return
   end
   
   [ hds, message ] = getmhs( mh.file_system );
   if isempty( message )
      message = writehd( fidphys, mh, hds )
   end
   fclose( fidphys );
   if ~isempty( message )
      res1 = message;
      return
   end
   pers.filelist{ pers.handlelist( fid )  }.mainheader = mh;
   return
   
case 'mh'
   % [ mh, message ] = ecatfile( 'mh', fid );
   res1 = [];
   res2 = '';
   fid = p1;
   if fid <= 0 | fid > length( pers.handlelist ) | fid <= 0
      res2 = 'Illegal handle';
      return;
   elseif pers.handlelist( fid ) == 0
      res2 = 'Illegal handle';
      return;
   end
   res1 = pers.filelist{ pers.handlelist( fid ) }.mainheader;
   return
   
case { 'open', 'create' }
% [ fid, message ]       = ecatfile( 'open', filename, file_system );
% [ fid, message ]       = ecatfile( 'open', filename, 'ecat7' );
% [ fid, message ]       = ecatfile( 'open', filename, 'ecat6.4' );
% [ fid, message ]       = ecatfile( 'open', filename, '' );
% [ fid, message ]       = ecatfile( 'create', filename, mh );
   filename = p1;
   if strcmp( filename, '' )
      fid = -1;
      message = 'Illegal filename';
      res1 = fid;
      res2 = message;
      return
   end
   
   if nargin < 3
      p2 = '';
   end
   
   for i = 1:length( pers.filelist )
      if ~isempty( pers.filelist{ i } )
         if strcmp( pers.filelist{ i }.filename, filename )
            if strcmp( action, 'create' )
               % file already open. Truncate the file.
               [ file, message ] = openecat( filename, 'create', p2 );
               if ~isempty( message )
                  fid = -1;
                  res1 = fid;
                  res2 = message;
                  return
               end
               pers.filelist{ i } = file;
            end
            % file already open, return a duplicate fid
            fid = firstno( [ pers.handlelist, 0 ] == 0 );
            pers.handlelist( fid ) = i;
            message = '';
            res1 = fid;
            res2 = message;
            return
         end
      end
   end
   
   [ file, message ] = openecat( filename, action, p2 );
   if ~isempty( message )
      fid = -1;
      res1 = fid;
      res2 = message;
      return
   end
   
   pos = length( pers.filelist ) + 1;
   for i = 1:length( pers.filelist )
      if isempty( pers.filelist{ i } ) 
         pos = i;
         break
      end
   end
   
   pers.filelist{ pos } = file;
   fid = firstno( [ pers.handlelist, 0 ] == 0 );
   pers.handlelist( fid ) = pos;
   message = '';
   res1 = fid;
   res2 = message;
   return
   
otherwise
   error( 'Illegal action' )
   
end % switch action

% ------------------------------------------------------------
function matnum = fpgdb2matnum( fpgdb )
frame = fpgdb( 1 );
plane = fpgdb( 2 );
gate = fpgdb( 3 );
data = fpgdb( 4 );
bed = fpgdb( 5 );
matnum = bitand( frame, 511 ) + ...
   bitshift( bitand( plane, 3 * 256 ), 9 - 8 ) + ...
   bitshift( bitand( data, 4 ), 11 - 2 ) + ...
   bitshift( bitand( bed, 15 ), 12 ) + ...
   bitshift( bitand( plane, 255 ), 16 ) + ...
   bitshift( bitand( gate, 63 ), 24 ) + ...
   bitshift( bitand( data, 3 ), 30 );

% ------------------------------------------------------------

function [ file, message ] = openecat( filename, create, file_system );
% function [ file, message ] = openecat( filename, 'open', file_system );
% function [ file, message ] = openecat( filename, 'create', mh );
if nargin < 2
   create = 'open';
end

message = '';
file = [];

% read mh, matdir, fpgdb

if strcmp( create, 'create' )
   % truncate or create
   mh = file_system;
   file_system = mh.file_system;
   switch mh.file_system
   case 'ecat7'
      fidphys = fopen( filename, 'w', 'ieee-be'  );
   case 'ecat6.4'
      fidphys = fopen( filename, 'w', 'vaxg'  );
   otherwise
      error( 'Wrong header' )
   end
   if fidphys == -1
      message = [ 'Can''t open file: ', filename ];
      return
   end
   [ hds, message ] = getmhs( mh.file_system );
   if isempty( message )
      message = writehd( fidphys, mh, hds );
   end
   if ~isempty( message )
      fclose( fidphys );
      return
   end
   
   fwrite( fidphys, 31, 'uint32' ); % 31 unused entries
   fwrite( fidphys,  2, 'uint32' ); % 1: next index block
   fwrite( fidphys,  2, 'uint32' ); % 2: previous index block ? 0 
   fwrite( fidphys,  0, 'uint32' ); % 3: number of used entries in the block?
   fwrite( fidphys, zeros( 4 * 31, 1 ), 'uint32' );
   if ftell( fidphys ) ~= 1024
      message = 'File write error';
   end
   fclose( fidphys );
end % create or truncate

[ file, message ] = readmh( filename, file_system );

if ~isempty( message )
   file = [];
   return
end

if strcmp( create, 'create' )
   fidphys = fopen( filename, 'r+', file.mhmachineformat  ); % 'w+' ?
else
   fidphys = fopen( filename, 'r', file.mhmachineformat  ); % 'w+' ?
end

if fidphys == -1
   file = [];
   message = [ 'Can''t open file: ', filename ];
   return
end

% read indexblocks

% first entry of an index block:
% 0: number of unused entries?
% 1: next index block
% 2: previous index block
% 3: number of used entries in the block?

% succeeding entries in an index block:
% 0: matnum
% 1: first data block (subheader)
% 2: last data block
% 3: status: 1:exists; 0: allocated, no data yet written, -1: matrix deleted

fseek( fidphys, 512, -1 );

blockno = 2;
matdir = [ ];
ok = 1;

file.dirblocks.data = [ ];
file.dirblocks.filepos = [ ];

while ok
   fseek( fidphys, 512 * ( blockno - 1), -1 );
   [ a1, count]= fread( fidphys, 128, 'uint32');
   if count < 128
      file = [];
      message = [ 'Can''t read index block: ', filename ];
      return
   end
   
   file.dirblocks.data = [ file.dirblocks.data, a1 ];
   file.dirblocks.filepos = [ file.dirblocks.filepos, blockno ];
   
   if a1(2) == 2
      ok = 0;
   elseif any( a1(2) == file.dirblocks.filepos( 2:end ) ) | a1( 2 ) <= 0
      % error if the directory block numbers are circular or not positive
      message = 'Illegal index block numbers';
      file = [];
      return;
   end
   
   blockno = a1(2);
   
end

% matdir(1:8, :)

%  fclose(fidphys);

sz = size( file.dirblocks.data );

direntries = [ ];
file.pos.dirblock = [ ];
file.pos.relpos = [ ];
for i = 1:sz( 2 )
   jj = file.dirblocks.data( 4, i );
   if jj < 0 | jj > 31
      message = 'Wrong index block';
      file = [];
      return
   end
   direntries = [ direntries, reshape( file.dirblocks.data( 5:4 + jj * 4, i ), 4, jj ) ];
   file.pos.dirblock = [ file.pos.dirblock, i + zeros( 1, jj ) ];
   file.pos.relpos = [ file.pos.relpos, 1:jj ];
end

fpgdb = file.dirblocks.data( 5:128, : );

fpgdb = [ zeros(7, sz(2)*31); reshape( fpgdb, 4, sz(2)*31 ) ];

m = direntries( 1, : );
frame   = rem( m, 512 ); m = floor( m / 512 );
hiplane = rem( m, 4 ); m = floor( m / 4 );
hidata  = rem( m, 2 ); m = floor( m / 2 );
bed     = rem( m, 16 ); m = floor( m / 16 );
plane   = rem( m, 256 ) + hiplane * 256; m = floor( m / 256 );
gate    = rem( m, 64 ); m = floor( m / 64 );
data    = rem( m, 4 ) + hidata * 4; m = floor( m / 4 );

file.pos.fpgdb = [ frame; plane; gate; data; bed ];
file.pos.blockstart = direntries( 2, : );
file.pos.blockend = direntries( 3, : );
file.pos.status = direntries( 4, : );

fclose( fidphys );

% ------------------------------------------------------
function [ ranges, message ] = getranges( file )

ranges = [ ];
message = '';

% check homogeneity
if length( file.pos.fpgdb( 1, : ) ) == 1
   m1 = file.pos.fpgdb;
   m2 = m1;
else
   m1 = min( file.pos.fpgdb' )';
   m2 = max( file.pos.fpgdb' )';
end
tot = m2 - m1 + 1;
f = 1;
for i = 1:5
   f( i, 1 ) = f( end ) * tot( i );
end

if f(5) ~= length( file.pos.status )
   message = 'Inhomogeneous fpgdb\n';
   return
end
testexists = zeros( f(5), 1 );
p = 1 + [ 1, f(1:4 )'] * ...
   ( file.pos.fpgdb - m1( :, ones( 1, size( file.pos.fpgdb, 2 ) ) ) );

testexists( p ) = 1;

if sum(testexists) ~= length( file.pos.status )
   message = 'Inhomogeneous fpgdb\n';
   return
end

ranges = [ m1, m2 ]';
if any( file.pos.status ~= 1 )
   message = 'Matrices with wrong status\n';
end

return

% ------------------------------------------------------

function [ hd, message ] = blankhd( hds )
message = '';
hd = [];
n = 0;
pos = [];

for i =1:length( hds.s )
   pos( end + 1 ) = n;
   switch hds.s(i).type
   case 'int16'
      n = n + 2 * hds.s(i).noelements;
      hd = setfield( hd, hds.s(i).fieldname, zeros( 1, hds.s(i).noelements ) );
   case { 'float32', 'mat1st', 'mat2nd', 'int32' }
      n = n + 4 * hds.s(i).noelements;
      hd = setfield( hd, hds.s(i).fieldname, zeros( 1, hds.s(i).noelements ) );
   case 'char'
      n = n + 1 * hds.s(i).noelements;
      hd = setfield( hd, hds.s(i).fieldname, '' );
   end % switch hds.s(i).type
end
if n ~= hds.header_size
   error( 'Blank header, wrong internal description' )
end

% make sure that ecat 7 mh contains a non-zero magic number
if strcmp( hds.file_system, 'ecat7' )
   hd.magic_number = 'MATRIX7';   
end

% ---------------------------------------------

function message = writehd( fidphys, hd, hds )
message = '';
initialpos = ftell( fidphys );

% make sure that ecat 7 mh contains a non-zero magic number
if strcmp( hd.file_system, 'ecat7' )
   if isempty( hd.magic_number )
      hd.magic_number = 'MATRIX7'; 
   elseif all( hd.magic_number == 0 )
      hd.magic_number = 'MATRIX7'; 
   end
end

ok = 1;
for i =1:length( hds.s )
   val = getfield( hd, hds.s(i).fieldname );
   if strcmp( hds.s(i).type, 'char' )
      if length( val(:) ) < hds.s(i).noelements
         val = [ val(:); char( ...
               zeros( hds.s(i).noelements - length( val(:) ), 1 ) ) ];
      end
   end
   if strcmp( hds.s(i).type, 'mat1st' )
      if length( val(:) ) ~= 12
         error( 'illegal header length\n' )
      end
   elseif strcmp( hds.s(i).type, 'mat2nd' )
      % the size has already been tested
   elseif length( val(:) ) ~= hds.s(i).noelements
      error( 'illegal header length\n' )
   end
   
   switch hds.s(i).type 
   case 'int16'
      ok = ok & length( val(:) ) == fwrite( fidphys, val, 'int16' );
   case 'int32'
      ok = ok & length( val(:) ) == fwrite( fidphys, val, 'int32' );
   case 'float32'
      ok = ok & length( val(:) ) == fwrite( fidphys, val, 'float32' );
   case 'mat1st'
      val = val( 1:3, 1:3 )';
      ok = ok & length( val(:) ) == fwrite( fidphys, val(:), 'float32' );
   case 'mat2nd'
      val = val(:,4);
      ok = ok & length( val(:) ) == fwrite( fidphys, val, 'float32' );
   case 'char'
      ok = ok & length( val(:) ) == fwrite( fidphys, val, 'uchar' );
   end % switch hds.s(i).type
end
if ~ok
   message = 'Can''t write header';
   return
end

if ftell( fidphys ) - initialpos ~= hds.header_size
   error( 'Writing header, wrong internal description' )
end

% ---------------------------------------------------
function [ file, message ] = readmh( filename, file_system )
% this function determines the header type and the machineformat 
% that is used for reading the main header and the index blocks.
% if the file is too short, it returns file = [], otherwise, it returns
%   file.mainheader = mh;
%   file.mhmachineformat = mhmachineformat;
%   file.filename = filename;
% the routine uses ad hoc testing to determine the format. This testing may
% need to be changed if the routine can not determine the format correctly   

message = '';
file = [];
mhmachineformat = 'ieee-be';
mh = [];

fidphys = fopen( filename, 'r', mhmachineformat );
if fidphys == -1
   message = [ 'Can''t open file: ', filename ];
   return
end

magic_number = fread( fidphys, 14, 'char' );
% check length of main header
[ mh, count ]= fread( fidphys, 256-7, 'uint16');
if count < 256-7
   message = 'Header read error';
   return
end

if ( all( magic_number == 0 ) & strcmp( file_system, '' ) ) | ...
      strcmp( file_system, 'ecat6.4' )
   fclose( fidphys );
   % ecat 6.4
   mhmachineformat = 'vaxg';
   fidphys = fopen( filename, 'r', mhmachineformat );
   if fidphys == -1
      return
   end
   hds = getmhs( 'ecat6.4' );
elseif strcmp( file_system, '' ) | strcmp( file_system, 'ecat7' ) 
   % ecat7
   fseek( fidphys, 0, -1 );
   hds = getmhs( 'ecat7' );
else
   message = 'Illegal file system';
   fclose( fidphys );
   return
end

[ mh, message ] = readhd( fidphys, hds );
if ~isempty( mh )
   file.mainheader = mh;
   file.mhmachineformat = mhmachineformat;
   file.filename = filename;
end

fclose( fidphys );

% -----------------------------------------------------

function [ hd, message ] = readhd( fidphys, hds )
message = '';
ok = 1;
hd = [];
initialpos = ftell( fidphys );

hd.file_system = hds.file_system;

for i =1:length( hds.s )
   switch hds.s(i).type
   case 'int16'
      [ val, count ] = fread( fidphys, hds.s(i).noelements, 'int16' );
      hd = setfield( hd, hds.s(i).fieldname, val' );
   case 'int32'
      [ val, count ] = fread( fidphys, hds.s(i).noelements, 'int32' );
      hd = setfield( hd, hds.s(i).fieldname, val' );
   case { 'float32', 'mat1st' }
      [ val, count ] = fread( fidphys, hds.s(i).noelements, 'float32' );
      hd = setfield( hd, hds.s(i).fieldname, val' );
   case 'mat2nd'
      [ val, count ] = fread( fidphys, hds.s(i).noelements, 'float32' );
      mt = getfield( hd, hds.s(i).fieldname );
      mt = [ reshape( mt, 3, 3 )', val ];
      hd = setfield( hd, hds.s(i).fieldname, mt );
   case 'char'
      [ val, count ] = fread( fidphys, hds.s(i).noelements, 'uchar' );
      val = char( val( val ~= 0 ) );
      hd = setfield( hd, hds.s(i).fieldname, val' );
   otherwise
      count = hds.s(i).noelements;
   end % switch hds.s(i).type
   ok = ok & count == hds.s(i).noelements;
   
end
if ~ok
   message = 'Header read error';
end
if ftell( fidphys ) - initialpos ~= hds.header_size
   error( 'Reading header, wrong internal description' )
end

function [ hds, message ] = getmhs( file_system )
message = '';

switch file_system
   
case 'ecat7'
   a = { 'magic_number',                   14, 'char'
      'original_file_name',             32, 'char'
      'sw_version',                      1, 'int16'
      'system_type',                     1, 'int16'
      'file_type',                       1, 'int16'
      'serial_number',                  10, 'char'
      'scan_start_time',                 1, 'int32'
      'isotope_name',                    8, 'char'
      'isotope_halflife',                1, 'float32'
      'radiopharmaceutical',            32, 'char'
      'gantry_tilt',                     1, 'float32'
      'gantry_rotation',                 1, 'float32'
      'bed_elevation',                   1, 'float32'
      'intrinsic_tilt',                  1, 'float32'
      'wobble_speed',                    1, 'int16'
      'transm_source_type',              1, 'int16'
      'distance_scanned',                1, 'float32'
      'transaxial_fov',                  1, 'float32'
      'angular_compression',             1, 'int16'
      'coin_samp_mode',                  1, 'int16'
      'axial_samp_mode',                 1, 'int16'
      'ecat_calibration_factor',         1, 'float32'
      'calibration_units',               1, 'int16'
      'calibration_units_label',         1, 'int16'
      'compression_code',                1, 'int16'
      'study_type',                     12, 'char'
      'patient_id',                     16, 'char'
      'patient_name',                   32, 'char'
      'patient_sex',                     1, 'char'
      'patient_dexterity',               1, 'char'
      'patient_age',                     1, 'float32'
      'patient_height',                  1, 'float32'
      'patient_weight',                  1, 'float32'
      'patient_birth_date',              1, 'int32'
      'physician_name',                 32, 'char'
      'operator_name',                  32, 'char'
      'study_description',              32, 'char'
      'acquisition_type',                1, 'int16'
      'patient_orientation',             1, 'int16'
      'facility_name',                  20, 'char'
      'num_planes',                      1, 'int16'
      'num_frames',                      1, 'int16'
      'num_gates',                       1, 'int16'
      'num_bed_pos',                     1, 'int16'
      'bed_position',                   16, 'float32'
      'plane_separation',                1, 'float32'
      'lwr_sctr_thres',                  1, 'int16'
      'lwr_true_thres',                  1, 'int16'
      'upr_true_thres',                  1, 'int16'
      'user_process_code',              10, 'char'
      'acquisition_mode',                1, 'int16'
      'bin_size',                        1, 'float32'
      'branching_fraction',              1, 'float32'
      'dose_start_time',                 1, 'int32'
      'dosage',                          1, 'float32'
      'well_counter_corr_factor',        1, 'float32'
      'data_units',                     32, 'char'
      'septa_state',                     1, 'int16'
      'fill',                            6, 'int16' };
   
   hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
      'type', a(:,3) );
   hds.header_size = 512;
   hds.file_system = 'ecat7';
   
case 'ecat6.4'
   a = { 'fill1',                          14, 'int16' % changed type
      'original_file_name',             20, 'char' % changed length
      'sw_version',                      1, 'int16'
      'data_type',                       1, 'int16' % ecat6.4
      'system_type',                     1, 'int16'
      'file_type',                       1, 'int16'
      'node_id',                        10, 'char' % serial_number
      ... %'scan_start_time', 1, 'int32'
      'scan_start_day',                  1, 'int16' % ecat6.4
      'scan_start_month',                1, 'int16' % ecat6.4
      'scan_start_year',                 1, 'int16' % ecat6.4
      'scan_start_hour',                 1, 'int16' % ecat6.4
      'scan_start_minute',               1, 'int16' % ecat6.4
      'scan_start_second',               1, 'int16' % ecat6.4
      'isotope_code',                    8, 'char' % isotope_name
      'isotope_halflife',                1, 'float32'
      'radiopharmaceutical',            32, 'char'
      'gantry_tilt',                     1, 'float32'
      'gantry_rotation',                 1, 'float32'
      'bed_elevation',                   1, 'float32'
      ... %'intrinsic_tilt', 1, 'float32'
      'rot_source_speed',                1, 'int16' % ecat6.4
      'wobble_speed',                    1, 'int16'
      'transm_source_type',              1, 'int16'
      ... %'distance_scanned', 1, 'float32'
      'axial_fov',                       1, 'float32'
      'transaxial_fov',                  1, 'float32'
      'transaxial_samp_mode',            1, 'int16' % ecat6.4
      ... %'angular_compression', 1, 'int16'
      'coin_samp_mode',                  1, 'int16'
      'axial_samp_mode',                 1, 'int16' % subtracted 1 in ecat7
      'calibration_factor',              1, 'float32' % ecat_calibration_factor
      'calibration_units',               1, 'int16'
      ... %'calibration_units_label', 1, 'int16'
      'compression_code',                1, 'int16'
      'study_name',                     12, 'char' % study_type
      'patient_id',                     16, 'char'
      'patient_name',                   32, 'char' % changed length
      'patient_sex',                     1, 'char'
      'patient_age',                    10, 'char' % changed type
      'patient_height',                 10, 'char' % changed type
      'patient_weight',                 10, 'char' % changed type
      'patient_dexterity',               1, 'char' % changed position
      ... %'patient_birth_date', 1, 'int32'
      'physician_name',                 32, 'char'
      'operator_name',                  32, 'char'
      'study_description',              32, 'char'
      'acquisition_type',                1, 'int16'
      'bed_type',                        1, 'int16' % ecat6.4
      'septa_type',                      1, 'int16'
      ... %'patient_orientation', 1, 'int16'
      'facility_name',                  20, 'char'
      'num_planes',                      1, 'int16'
      'num_frames',                      1, 'int16'
      'num_gates',                       1, 'int16'
      'num_bed_pos',                     1, 'int16'
      'bed_position',                   16, 'float32'
      'plane_separation',                1, 'float32'
      'lwr_sctr_thres',                  1, 'int16'
      'lwr_true_thres',                  1, 'int16'
      'upr_true_thres',                  1, 'int16'
      'collimator',                      1, 'float32' % ecat6.4
      'user_process_code',              10, 'char'
      ... %'acquisition_mode', 1, 'int16'
      'fill2',                          20, 'int16' };
      %'bin_size', 1, 'float32'
      %'branching_fraction', 1, 'float32'
      %'dose_start_time', 1, 'int32'
      %'dosage', 1, 'float32'
      %'well_counter_corr_factor', 1, 'float32'
      %'data_units', 32, 'char'
      %'septa_state', 1, 'int16'
      %'fill', 6, 'int16' };
   
   hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
      'type', a(:,3) );
   hds.header_size = 512;
   hds.file_system = 'ecat6.4';
otherwise
   message = 'Illegal file system';
   hds = [];

end % switch file_system

function [ hds, message ] = getshs( file_system, file_type );
message = '';
hds = [];

switch file_system
case 'ecat7'
   
   switch file_type
      
      % implemented:
      % 01=Sinogram, 03=Attenuation Correction, 04=Normalization, 
      % 05=Polar Map, 07=Volume 16, 11=3D Sinogram 16, 13=3D Normalization
      
      % not implemented:
      % 00=unknown, 02=Image-16, 06=Volume 8, 08=Projection 8, 
      % 09=Projection 16, 10=Image 8, 12=3D Sinogram 8, 14=3D Sinogram Fit
      
   case 1
      % sinogram
      a = { ...
            'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'num_r_elements',                  1, 'int16'
         'num_angles',                      1, 'int16'
         'corrections_applied',             1, 'int16'
         'num_z_elements',                  1, 'int16'
         'ring_difference',                 1, 'int16'
         'xyz_resolution',                  3, 'float32'
         'w_resolution',                    1, 'float32'
         'fill1',                           6, 'int16'
         'gate_duration',                   1, 'int32'
         'r_wave_offset',                   1, 'int32'
         'num_accepted_beats',              1, 'int32'
         'scale_factor',                    1, 'float32'
         'scan_min',                        1, 'int16'
         'scan_max',                        1, 'int16'
         'prompts',                         1, 'int32'
         'delayed',                         1, 'int32'
         'multiples',                       1, 'int32'
         'net_trues',                       1, 'int32'
         'cor_singles',                    16, 'float32'
         'uncor_singles',                  16, 'float32'
         'tot_avg_cor',                     1, 'float32'
         'tot_avg_uncor',                   1, 'float32'
         'total_coin_rate',                 1, 'int32'
         'frame_start_time',                1, 'int32'
         'frame_duration',                  1, 'int32'
         'deadtime_correction_factor',      1, 'float32'
         'physical_planes',                 8, 'int16'
         'fill2',                          83, 'int16'
         'fill3',                          50, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
   case 3
      % attenuation
      a = { ...
            'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'attenuation_type',                1, 'int16'
         'num_r_elements',                  1, 'int16'
         'num_angles',                      1, 'int16'
         'num_z_elements',                  1, 'int16'
         'ring_difference',                 1, 'int16'
         'x_resolution',                    1, 'float32'
         'y_resolution',                    1, 'float32'
         'z_resolution',                    1, 'float32'
         'w_resolution',                    1, 'float32'
         'scale_factor',                    1, 'float32'
         'x_offset',                        1, 'float32'
         'y_offset',                        1, 'float32'
         'x_radius',                        1, 'float32'
         'y_radius',                        1, 'float32'
         'tilt_angle',                      1, 'float32'
         'attenuation_coeff',               1, 'float32'
         'attenuation_min',                 1, 'float32'
         'attenuation_max',                 1, 'float32'
         'skull_thickness',                 1, 'float32'
         'num_additional_atten_coeff',      1, 'int16'
         'additional_atten_coeff',          8, 'float32'
         'edge_finding_threshold',          1, 'float32'
         'storage_order',                   1, 'int16'
         'span',                            1, 'int16'
         'z_elements',                     64, 'int16'
         'fill1',                          86, 'int16'
         'fill2',                          50, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
   case 4
      % 6.5 normalization
      a = { ...
            'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'num_r_elements',                  1, 'int16'
         'num_angles',                      1, 'int16'
         'num_z_elements',                  1, 'int16'
         'ring_difference',                 1, 'int16'
         'scale_factor',                    1, 'float32'
         'norm_min',                        1, 'float32'
         'norm_max',                        1, 'float32'
         'fov_source_width',                1, 'float32'
         'norm_quality_factor',             1, 'float32'
         'norm_quality_factor_code',        1, 'int16'
         'storage_order',                   1, 'int16'
         'span',                            1, 'int16'
         'z_elements',                     64, 'int16'
         'fill1',                         123, 'int16'
         'fill2',                          50, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
      % polar map
   case 5
      a = { ...
            'data_type',                       1, 'int16'
         'polar_map_type',                  1, 'int16'
         'num_rings',                       1, 'int16'
         'sectors_per_ring',               32, 'int16'
         'ring_position',                  32, 'float32'
         'ring_angle',                     32, 'int16'
         'start_angle',                     1, 'int16'
         'long_axis_left',                  3, 'int16'
         'long_axis_right',                 3, 'int16'
         'position_data',                   1, 'int16'
         'image_min',                       1, 'int16'
         'image_max',                       1, 'int16'
         'scale_factor',                    1, 'float32'
         'pixel_size',                      1, 'float32'
         'frame_duration',                  1, 'int32'
         'frame_start_time',                1, 'int32'
         'processing_code',                 1, 'int16'
         'quant_units',                     1, 'int16'
         'annotation',                     40, 'char'
         'gate_duration',                   1, 'int32'
         'r_wave_offset',                   1, 'int32'
         'num_accepted_beats',              1, 'int32'
         'polar_map_protocol',             20, 'char'
         'database_name',                  30, 'char'
         'fill1',                          27, 'int16'
         'fill2',                          27, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
   case 7
      % image
      a = { ...
            'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'xyz_dimension',                   3, 'int16'
         'xyz_offset',                      3, 'float32'
         'recon_zoom',                      1, 'float32'
         'scale_factor',                    1, 'float32'
         'image_min',                       1, 'int16'
         'image_max',                       1, 'int16'
         'xyz_pixel_size',                  3, 'float32'
         'frame_duration',                  1, 'int32'
         'frame_start_time',                1, 'int32'
         'filter_code',                     1, 'int16'
         'xyz_resolution',                  3, 'float32'
         'num_r_elements',                  1, 'float32'
         'num_angles',                      1, 'float32'
         'z_rotation_angle',                1, 'float32'
         'decay_corr_fctr',                 1, 'float32'
         'processing_code',                 1, 'int32'
         'gate_duration',                   1, 'int32'
         'r_wave_offset',                   1, 'int32'
         'num_accepted_beats',              1, 'int32'
         'filter_cutoff_frequency',         1, 'float32'
         'filter_resolution',               1, 'float32'
         'filter_ramp_slope',               1, 'float32'
         'filter_order',                    1, 'int16'
         'filter_scatter_fraction',         1, 'float32'
         'filter_scatter_slope',            1, 'float32'
         'annotation',                     40, 'char'
         'mt',                              9, 'mat1st'
         'rfilter_cutoff',                  1, 'float32'
         'rfilter_resolution',              1, 'float32'
         'rfilter_code',                    1, 'int16'
         'rfilter_order',                   1, 'int16'
         'zfilter_cutoff',                  1, 'float32'
         'zfilter_resolution',              1, 'float32'
         'zfilter_code',                    1, 'int16'
         'zfilter_order',                   1, 'int16'
         'mt',                              3, 'mat2nd'
         'scatter_type',                    1, 'int16'
         'recon_type',                      1, 'int16'
         'recon_views',                     1, 'int16'
         'fill',                          136, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
   case { 11, 14 }
      % 3D sinogram
      a = { ...
            'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'num_r_elements',                  1, 'int16'
         'num_angles',                      1, 'int16'
         'corrections_applied',             1, 'int16'
         'num_z_elements',                 64, 'int16'
         'ring_difference',                 1, 'int16'
         'storage_order',                   1, 'int16'
         'axial_compression',               1, 'int16'
         'x_resolution',                    1, 'float32'
         'v_resolution',                    1, 'float32'
         'z_resolution',                    1, 'float32'
         'w_resolution',                    1, 'float32'
         'fill1',                           6, 'int16'
         'gate_duration',                   1, 'int32'
         'r_wave_offset',                   1, 'int32'
         'num_accepted_beats',              1, 'int32'
         'scale_factor',                    1, 'float32'
         'scan_min',                        1, 'int16'
         'scan_max',                        1, 'int16'
         'prompts',                         1, 'int32'
         'delayed',                         1, 'int32'
         'multiples',                       1, 'int32'
         'net_trues',                       1, 'int32'
         'tot_avg_cor',                     1, 'float32'
         'tot_avg_uncor',                   1, 'float32'
         'total_coin_rate',                 1, 'int32'
         'frame_start_time',                1, 'int32'
         'frame_duration',                  1, 'int32'
         'deadtime_correction_factor',      1, 'float32'
         'fill2',                          90, 'int16'
         'fill3',                          50, 'int16'
         'uncor_singles',                 128, 'float32' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 1024;
      
   case 13
      % 3D normalization
      a = { ...
            'data_type',                       1, 'int16'
         'num_r_elements',                  1, 'int16'
         'num_transaxial_crystals',         1, 'int16'
         'num_crystal_rings',               1, 'int16'
         'crystals_per_ring',               1, 'int16'
         'num_geo_corr_planes',             1, 'int16'
         'uld',                             1, 'int16'
         'lld',                             1, 'int16'
         'scatter_energy',                  1, 'int16'
         'norm_quality_factor',             1, 'float32'
         'norm_quality_factor_code',        1, 'int16'
         'ring_dtcor1',                    32, 'float32'
         'ring_dtcor2',                    32, 'float32'
         'crystal_dtcor',                   8, 'float32'
         'span',                            1, 'int16'
         'max_ring_diff',                   1, 'int16'
         'fill1',                          48, 'int16'
         'fill2',                          50, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
   otherwise
      message = 'Illegal file type';
      
   end % switch file_type
   
   hds.file_system = 'ecat7sh';
   
case 'ecat6.4'
   
   switch file_type
      % implemented:
      % 01=.scn, 02=.img, 03=.atn, 04=.nrm 
      
   case 1 % .scn
      
      a = { ...
            'fill1',                          63, 'int16'
         'data_type',                       1, 'int16'
         'fill2',                           2, 'int16'
         'dimension_1',                     1, 'int16'
         'dimension_2',                     1, 'int16'
         'smoothing',                       1, 'int16'
         'processing_code',                 1, 'int16'
         'fill3',                           3, 'int16'
         'sample_distance',                 1, 'float32'
         'fill4',                           8, 'int16'
         'isotope_halflife',                1, 'float32'
         'frame_duration_sec',              1, 'int16'
         'gate_duration',                   1, 'int32'
         'r_wave_offset',                   1, 'int32'
         'fill5',                           1, 'int16'
         'scale_factor',                    1, 'float32'
         'fill6',                           3, 'int16'
         'scan_min',                        1, 'int16'
         'scan_max',                        1, 'int16'
         'prompts',                         1, 'int32'
         'delayed',                         1, 'int32'
         'multiples',                       1, 'int32'
         'net_trues',                       1, 'int32'
         'fill7',                          52, 'int16'
         'cor_singles',                    16, 'float32'
         'uncor_singles',                  16, 'float32'
         'tot_avg_cor',                     1, 'float32'
         'tot_avg_uncor',                   1, 'float32'
         'total_coin_rate',                 1, 'int32'
         'frame_start_time',                1, 'int32'
         'frame_duration',                  1, 'int32'
         'loss_correction_fctr'             1, 'float32'
         'fill8',                          22, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
      
   case 2 % .img
      
      a = { ...
            'fill1',                          63, 'int16'
         'data_type',                       1, 'int16'
         'num_dimensions',                  1, 'int16'
         'fill2',                           1, 'int16'
         'dimension_1',                     1, 'int16'
         'dimension_2',                     1, 'int16'
         'fill3',                          12, 'int16'
         'x_origin',                        1, 'float32'
         'y_origin',                        1, 'float32'
         'recon_scale',                     1, 'float32'
         'quant_scale',                     1, 'float32'
         'image_min',                       1, 'int16'
         'image_max',                       1, 'int16'
         'fill4',                           2, 'int16'
         'pixel_size',                      1, 'float32'
         'slice_width',                     1, 'float32'
         'frame_duration',                  1, 'int32'
         'frame_start_time',                1, 'int32'
         'slice_location',                  1, 'int16'
         'recon_start_hour',                1, 'int16'
         'recon_start_min',                 1, 'int16'
         'recon_start_sec',                 1, 'int16'
         'recon_duration',                  1, 'int32'
         'fill5',                          12, 'int16'
         'filter_code',                     1, 'int16'
         'scan_matrix_num',                 1, 'int32'
         'norm_matrix_num',                 1, 'int32'
         'atten_cor_mat_num',               1, 'int32'
         'fill6',                          23, 'int16'
         'image_rotation',                  1, 'float32'
         'plane_eff_corr_fctr',             1, 'float32'
         'decay_corr_fctr',                 1, 'float32'
         'loss_corr_fctr',                  1, 'float32'
         'fill7',                          32, 'int16'
         'processing_code',                 1, 'int16'
         'fill8',                           1, 'int16'
         'quant_units',                     1, 'int16'
         'recon_start_day',                 1, 'int16'
         'recon_start_month',               1, 'int16'
         'recon_start_year',                1, 'int16'
         'ecat_calibration_fctr',           1, 'float32'
         'well_counter_cal_fctr',           1, 'float32'
         'filter_params',                   6, 'float32'
         'annotation',                     40, 'char'
         'fill9',                          26, 'int16' };
         
         hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
            'type', a(:,3) );
         hds.header_size = 512;
      
   case 3 % .atn
      
      a = { ...
            'fill1',                          63, 'int16'
         'data_type',                       1, 'int16'
         'attenuation_type',                1, 'int16'
         'fill2',                           1, 'int16'
         'dimension_1',                     1, 'int16'
         'dimension_2',                     1, 'int16'
         'fill3',                          23, 'int16'
         'scale_factor',                    1, 'float32'
         'x_origin',                        1, 'float32'
         'y_origin',                        1, 'float32'
         'x_radius',                        1, 'float32'
         'y_radius',                        1, 'float32'
         'tilt_angle',                      1, 'float32'
         'attenuation_coeff',               1, 'float32'
         'sample_distance',                 1, 'float32'
         'fill4',                         149, 'int16' };
         
         hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
            'type', a(:,3) );
         hds.header_size = 512;
      
   case 4 % .nrm
      
      a = { ...
            'fill1',                          63, 'int16'
         'data_type',                       1, 'int16'
         'fill2',                           2, 'int16'
         'dimension_1',                     1, 'int16'
         'dimension_2',                     1, 'int16'
         'fill3',                          23, 'int16'
         'scale_factor',                    1, 'float32'
         'norm_hour',                       1, 'int16'
         'norm_minute',                     1, 'int16'
         'norm_second',                     1, 'int16'
         'norm_day',                        1, 'int16'
         'norm_month',                      1, 'int16'
         'norm_year',                       1, 'int16'
         'fov_source_width',                1, 'float32'
         'fill4',                         155, 'int16' };
      
      hds.s = struct( 'fieldname', a(:,1), 'noelements', a(:,2), ...
         'type', a(:,3) );
      hds.header_size = 512;
   otherwise
      message = 'Illegal file type';
      
   end % switch file_type
   
   hds.file_system = 'ecat6.4sh';
   
otherwise
   message = 'Illegal file system';
   
end % switch file_system

function [ datadescr, message ] = getdatadescr( file_system, file_type, data_type );

datadescr = [];
message = '';

switch file_system
   
case 'ecat6.4'
   switch file_type
      
   case { 1, 2, 3, 4 }
      % 1 = .scn, 2 = .img, 3 = .atn, 4 = .nrm
      switch data_type
      case 2 % VAX_Ix2
         datadescr.ftype = 'int16';
         datadescr.size = 2;
      case 4 % vax float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      case 5 % ieee float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      otherwise
         message = 'Data type not implemented';
      end
   otherwise 
      message = 'File type not implemented';
   end %switch file_type
   
case 'ecat7'
   
   % implemented:
   % 01=Sinogram, 03=Attenuation Correction, 04=Normalization, 
   % 05=Polar Map, 07=Volume 16, 11=3D Sinogram 16, 13=3D Normalization
   
   % not implemented:
   % 00=unknown, 02=Image-16, 06=Volume 8, 08=Projection 8, 
   % 09=Projection 16, 10=Image 8, 12=3D Sinogram 8, 14=3D Sinogram Fit
   switch file_type
   case { 11, 14 }
      switch data_type
      case 1 % byte
         datadescr.ftype = 'int8';
         datadescr.size = 1;
      case { 2, 6 } % Sun short, sun_in
         datadescr.ftype = 'int16';
         datadescr.size = 2;
      case 5 % ieee float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      otherwise
         message = 'Data type not implemented';
      end   
   case 13
      switch data_type
      case 1 % ieee-float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      case 5 % sun-float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      otherwise
         message = 'Data type not implemented';
      end   
   case { 1, 3, 4, 5, 7 }
      % (0=Unkonwn Matrix Data Type, 1=Byte Data, 2=VAX_Ix2, 3=VAX_Ix4, 
      %  4=VAX_Rx4, 5=IEEE Float, 6=Sun short, 7=Sun long) 
      % DTYPE_BYTES, _I2, _I4, _VAXR4, _SUNFL, _SUNIN 
      switch data_type    
      case 1 % byte
         datadescr.ftype = 'int8';
         datadescr.size = 1;
      case 6 % Sun short
         datadescr.ftype = 'int16';
         datadescr.size = 2;
      case 5 % ieee float
         datadescr.ftype = 'float32';
         datadescr.size = 4;
      otherwise
         message = 'Data type not implemented';
      end
   otherwise 
      error( 'File type not implemented' );
   end %switch file_type
   
otherwise
   message = 'File system not implemented';
   
end % switch file_system
% ------------------------------------------------------------

function [ hd, offsetfromheader, writefill, numdatablocks, message ] = ...
      set_xyz_dimension( hd, datadescr, xyz_dimension, selsegment  );
% implemented:
% 01=Sinogram, 03=Attenuation Correction, 04=Normalization,
% 05=Polar Map, 07=Volume 16, 11=3D Sinogram 16, 13=3D Normalization

message = '';
offsetfromheader = 0;
writefill = 0;
numdatablocks = 0;

switch hd.mh.file_system
case 'ecat6.4'
   
   switch hd.mh.file_type
   case { 01, 02, 03, 04 } % .scn, .img, .atn, .nrm
      hd.sh.dimension_1 = xyz_dimension( 1 );
      hd.sh.dimension_2 = xyz_dimension( 2 );
      numdatablocks = floor( ...
         ( prod( xyz_dimension ) * datadescr.size + 511 ) / 512 );
      writefill = ( 511 - rem( prod( xyz_dimension ) * datadescr.size + 511, 512 ) );
      
   otherwise
      message = 'Writing of the file type not implemented';
      hd = [];
      
   end
   
case 'ecat7'
   
   switch hd.mh.file_type
   case 07 % Volume 16
      hd.sh.xyz_dimension = xyz_dimension;
      numdatablocks = floor( ...
         ( prod( xyz_dimension ) * datadescr.size + 511 ) / 512 );
      writefill = ( 511 - rem( prod( xyz_dimension ) * datadescr.size + 511, 512 ) );
      
   case { 11, 14 } % 3D Sinogram 16, FORE sinogram float
      if hd.sh.storage_order == 0
         hd.sh.num_r_elements = xyz_dimension( 1 );
         hd.sh.num_z_elements( selsegment ) = xyz_dimension( 2 );
         hd.sh.num_angles = xyz_dimension( 3 );
      else
         hd.sh.num_r_elements = xyz_dimension( 1 );
         hd.sh.num_angles = xyz_dimension( 2 );
         hd.sh.num_z_elements( selsegment ) = xyz_dimension( 3 );
      end
      
      nobytes = hd.sh.num_r_elements * hd.sh.num_angles * datadescr.size * ...
         hd.sh.num_z_elements;
      %   numdatablocksarr = floor( ( nobytes + 511 ) / 512 );
      %   offsetfromheader = sum( numdatablocksarr( 1:selsegment-1 ) ) * 512;
      offsetfromheader = sum( nobytes( 1:selsegment-1 ) );
      
      %   numdatablocks = sum( numdatablocksarr );
      numdatablocks = floor( ( sum( nobytes ) + 511 ) / 512 );
      
      if sum( hd.sh.num_z_elements( selsegment + 1 : end ) ) == 0
         writefill = ( 511 - rem( prod( xyz_dimension ) * datadescr.size + offsetfromheader + 511, 512 ) );
      else
         writefill = 0;
      end
      
   case 3 % Attenuation
      if hd.sh.storage_order == 0
         hd.sh.num_r_elements = xyz_dimension( 1 );
         hd.sh.z_elements( selsegment ) = xyz_dimension( 2 );
         hd.sh.num_angles = xyz_dimension( 3 );
      else
         hd.sh.num_r_elements = xyz_dimension( 1 );
         hd.sh.num_angles = xyz_dimension( 2 );
         hd.sh.z_elements( selsegment ) = xyz_dimension( 3 );
      end
      if selsegment == 1
         hd.sh.num_z_elements = hd.sh.z_elements( 1 );
      end
      
      nobytes = hd.sh.num_r_elements * hd.sh.num_angles * datadescr.size * ...
         hd.sh.z_elements;
      %   numdatablocksarr = floor( ( nobytes + 511 ) / 512 );
      %   offsetfromheader = sum( numdatablocksarr( 1:selsegment-1 ) ) * 512;
      offsetfromheader = sum( nobytes( 1:selsegment-1 ) );
      
      %   numdatablocks = sum( numdatablocksarr );
      numdatablocks = floor( ( sum( nobytes ) + 511 ) / 512 );
      
      if sum( hd.sh.num_z_elements( selsegment + 1 : end ) ) == 0
         writefill = ( 511 - rem( prod( xyz_dimension ) * datadescr.size + offsetfromheader + 511, 512 ) );
      else
         writefill = 0;
      end
      
   otherwise
      message = 'Writing of the file type not implemented';
      hd = [];
      
   end
   
otherwise
   message = 'File system not implemented';
   hd = [];
   
end % switch hd.mh.file_system

% ------------------------------------------------------------
function [ xyz_dimension, offsetfromheader, message ] = get_xyz_dimension( hd, datadescr, selsegment )
% implemented:
% 01=Sinogram, 03=Attenuation Correction, 04=Normalization, 
% 05=Polar Map, 07=Volume 16, 11=3D Sinogram 16, 13=3D Normalization

message = '';
xyz_dimension = [];
offsetfromheader = 0;

switch hd.mh.file_system
case 'ecat7'
   
   switch hd.mh.file_type
   case 07 % Volume 16
      xyz_dimension = hd.sh.xyz_dimension;
      
   case { 11, 14 } % 3D Sinogram 16, FORE sinogram float
      if hd.sh.storage_order == 0
         xyz_dimension = [ hd.sh.num_r_elements, ...
               hd.sh.num_z_elements( selsegment ), hd.sh.num_angles ];
      else
         xyz_dimension = [ hd.sh.num_r_elements, hd.sh.num_angles, ...
               hd.sh.num_z_elements( selsegment ) ];
      end
      nobytes = hd.sh.num_r_elements * hd.sh.num_angles * datadescr.size * ...
         hd.sh.num_z_elements;
      %   numdatablocksarr = floor( ( nobytes + 511 ) / 512 );
      %   offsetfromheader = sum( numdatablocksarr( 1:selsegment-1 ) ) * 512;
      offsetfromheader = sum( nobytes( 1:selsegment-1 ) );
      
   case 3 % Attenuation
      if hd.sh.storage_order == 0
         xyz_dimension = [ hd.sh.num_r_elements, ...
               hd.sh.z_elements( selsegment ), hd.sh.num_angles ];
      else
         xyz_dimension = [ hd.sh.num_r_elements, hd.sh.num_angles, ...
               hd.sh.z_elements( selsegment ) ];
      end
      nobytes = hd.sh.num_r_elements * hd.sh.num_angles * datadescr.size * ...
         hd.sh.z_elements;
      %   numdatablocksarr = floor( ( nobytes + 511 ) / 512 );
      %   offsetfromheader = sum( numdatablocksarr( 1:selsegment-1 ) ) * 512;
      offsetfromheader = sum( nobytes( 1:selsegment-1 ) );
      
   otherwise
      message = 'Reading of the file type not implemented';
      
   end
   
case 'ecat6.4'
   
   switch hd.mh.file_type
   case { 01, 02, 03, 04 } % .scn, .img, .atn, .nrm
      xyz_dimension = [ hd.sh.dimension_1, hd.sh.dimension_2, 1 ];
      
   otherwise
      message = 'Reading of the file type not implemented';
      
   end
   
otherwise
   message = 'File system not implemented';
   
end % switch hd.mh.file_system

% ------------------------------------------------------------

function n = firstno( sel )
n = [ find( sel(:) ); 0 ];
n = n( 1 );

% ------------------------------------------------------------