function test_inq_att ( ncfile )
% TEST_INQ_ATT:  
%
% The matlab API is
%
%     [datatype, attlen, status] = mexnc ( 'inq_att', ncid, varid, attname );
%    
%
% Test 1:  Normal retrieval.
% Test 2:  Invalid ncid 
% Test 3:  Invalid varid 
% Test 4:  Invalid name. 
% Test 5:  ncid = []
% Test 6:  varid = [] 
% Test 7:  name = []
% Test 8:  non numeric ncid
% Test 9:  non numeric varid
% Test 10: non character attribute name


error_condition = 0;

[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end


[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', 20 );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end

[varid, status] = mexnc ( 'def_var', ncid, 'x', nc_double, 1, xdimid );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end

input_data = [3.14159 0];
status = mexnc ( 'put_att_double', ncid, varid, 'test_double', nc_double, 2, input_data );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end
[status] = mexnc ( 'enddef', ncid );

status = mexnc ( 'sync', ncid );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end


% Test 1
testid = 'Test 1';
[datatype, att_length, status] = mexnc ( 'inq_att', ncid, varid, 'test_double' );
if ( status ~= 0 )
	err_msg = sprintf ( '%s:  %s:  ''%s''\n', mfilename, ncerr, mexnc ( 'strerror', status ) );
	error ( err_msg );
end

if datatype ~= nc_double
	msg = sprintf ( '%s:  datatype did not match.\n', mfilename );
	error ( msg );
end

if att_length ~= 2
	msg = sprintf ( '%s:  attribute length did not match.\n', mfilename );
	error ( msg );
end



% Test 2:  Invalid ncid 
testid = 'Test 2';
[value,status] = mexnc ( 'INQ_ATT', -2000, varid, 'test_double' );
if ( status == 0 )
	err_msg = sprintf ( '%s:  %s:  ''%s''\n', mfilename, ncerr, mexnc ( 'strerror', status ) );
	error ( err_msg );
end

	
% Test 3:  Invalid varid 
testid = 'Test 3';
[value,status] = mexnc ( 'INQ_ATT', ncid,  -2000, 'test_double' );
if ( status == 0 )
	err_msg = sprintf ( '%s:  %s:  ''%s''\n', mfilename, ncerr, mexnc ( 'strerror', status ) );
	error ( err_msg );
end

	
% Test 4:  Invalid name. 
testid = 'Test 4';
[value,status] = mexnc ( 'INQ_ATT', ncid,  varid, 'test_double2' );
if ( status == 0 )
	err_msg = sprintf ( '%s:  %s:  ''%s''\n', mfilename, ncerr, mexnc ( 'strerror', status ) );
	error ( err_msg );
end

	
% Test 5:  ncid = []
testid = 'Test 5';
try
	[value,status] = mexnc ( 'INQ_ATT', [], varid, 'test_double' );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


% Test 6:  varid = [] 
testid = 'Test 6';
try
	[value,status] = mexnc ( 'INQ_ATT', ncid, [], 'test_double' );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


% Test 7:  name = []
testid = 'Test 7';
try
	[value,status] = mexnc ( 'INQ_ATT', ncid, varid, [] );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


% Test 8:  non numeric ncid
testid = 'Test 8';
try
	[value,status] = mexnc ( 'INQ_ATT', 'ncid', varid, 'test_double' );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


% Test 9:  non numeric varid
testid = 'Test 9';
try
	[value,status] = mexnc ( 'INQ_ATT', ncid, 'varid', 'test_double' );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


% Test 10: non character attribute name
testid = 'Test 10';
try
	[value,status] = mexnc ( 'INQ_ATT', ncid, varid, 0 );
	error_condition = 1;
end
if error_condition
	err_msg = sprintf ( '%s:  %s:  succeeded when it should have failed\n', mfilename, testid);
	error ( err_msg );
end


mexnc ( 'close', ncid );

fprintf ( 1, 'INQ_ATT succeeded.\n' );
return
