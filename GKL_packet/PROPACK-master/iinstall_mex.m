function install_mex( TESTING ) 
% run this code to compile mex files for PROPACK and SVT 
%       Stephen Becker, 2/21/09 srbecker@caltech.edu 
 
% For PROPACK: 
% Had to rewrite some of the fortran code into C, because 
% compiling fortran on Windows is a major headache. 
 
global NO_FORTRAN 
if isempty(NO_FORTRAN), NO_FORTRAN = false; end 
 
 
if nargin < 1 
TESTING = false;  % controls whether it compiles with the verbose ('-v') flag 
end 
 
if NO_FORTRAN 
    TESTING = true; 
end 
 
% Determine if the version is newer than 2006a, 
% in which case mex should be given the -largeArrayDims flag 
v = version('-release'); 
yr = str2num(v(1:4)); 
if yr > 2006 || ( yr == 2006 && v(5) == 'b' ) 
    R2006a = false; 
else 
    R2006a = true; 
end 
 
 
fprintf('\n'); 
in = input('This will compile mex files.  Continue? [y]/n  ','s'); 
if ~isempty(in) && (strcmp(in,'n') || strcmp(in,'no') ) 
    return; 
end 
 
 
c = computer; 
if ispc 
    cc = get_compiler_config(); 
 
    if strcmpi(cc,'microsoft') 
        cc = 'microsoft'; 
    else 
        cc = 'lcc';     % Matlab's bundled compiler 
    end 
    if strfind(c,'64') 
        disp('If you have 64-bit windows, you might need to modify this script yourself'); 
        libpath = fullfile(matlabroot,'extern','lib','win64',cc); 
    else 
        libpath = fullfile(matlabroot,'extern','lib','win32',cc); 
    end 
 
    % look for necessary files (i.e. make sure we're in right directory) 
    if ~exist(fullfile(pwd,'dbdqr.c'),'file') || ~exist(fullfile(pwd,'bdsqr_mex.c'),'file') 
        disp('Please run this from the same directory as dbdqr.c and bdsqr_mex.c'); 
        return; 
    end 
     
    % look for pre-built executables 
%     mexext('all');  c = computer  
    ext = mexext; 
    filef = @(name) fullfile(pwd,[name,'.',ext]); 
    INSTALL_PROPACK = true; 
    if exist(filef('bdsqr'),'file') 
        disp('Found precompiled PROPACK executables for your platform'); 
        in = input('Recompile anyway? (old executable will be copied to the "Old" subdirectory) [n]/y  ','s'); 
        if ~isempty(in) && ( strcmp(in,'y') || strcmp(in,'yes') ) 
            disp('Backing up old executables ...'); 
            if ~exist('Old','file'), mkdir('Old'); end 
            if ~exist(['Old',filesep,'bdsqr.',ext],'file') 
                copyfile(filef('bdsqr'),'Old'); 
            end 
            disp('Proceeding to compile...'); 
        else 
            in = input('Would you like to install the SVT executables?  [y]/n  ','s'); 
            if ~isempty(in) && ( strcmp(in,'n') || strcmp(in,'no') ) 
                return; 
            end 
            INSTALL_PROPACK = false; 
        end 
    end 
     
    % build a temporary file 
%     rspfile = [tempname '.rsp'];  % if this fails, try with  tempname(pwd) 
    rspfile = 'temp_file.rsp'; 
    [Frsp, errmsg] = fopen(rspfile, 'wt'); 
    if ~isempty(errmsg)  
%         rspfile = [tempname(pwd) '.rsp'];  
%         [Frsp, errmsg] = fopen(rspfile, 'wt'); 
%         if ~isempty(errmsg)  
            disp('Error: can''t open temp file');  return; 
%         end 
    end 
    fprintf(Frsp,'-L''%s''',libpath); 
    fclose(Frsp); 
     
    if INSTALL_PROPACK 
        % try to compile, using Matlab's default BLAS ("mwblas") 
        % do NOT compile reorth.c, because it is buggy! 
        if TESTING 
            mex -v -O -DWINDOWS @temp_file.rsp  -lmwlapack -lmwblas dbdqr.c bdsqr_mex.c -output bdsqr 
%             mex -v -O -DWINDOWS @temp_file.rsp  -lmwlapack -lmwblas reorth.c reorth_mex.c -output reorth 
        else 
            mex -O -DWINDOWS @temp_file.rsp  -lmwlapack -lmwblas dbdqr.c bdsqr_mex.c -output bdsqr 
%             mex -O -DWINDOWS @temp_file.rsp  -lmwlapack -lmwblas reorth.c reorth_mex.c -output reorth 
        end 
 
        % see if it worked 
        disp('Looks like compilation of PROPACK mex file worked.  Now, testing the installation'); 
        a = randn(50,1); b = randn(50,1); 
        [s,bnd]=bdsqr(a,b); 
        disp('No error messages!  Finished installing PROPACK.'); 
    end 
     
    % This is what we don't want to see: PROPACK:NotUsingMex','Using slow 
    % matlab code for bdsqr.' 
     
    % ---------- for SVT ------------ 
    disp('Now compiling SVT mex files'); 
    if exist(filef('XonOmega'),'file') || exist(filef('updateSparse'),'file') 
        disp('Found precompiled SVT executables for your platform'); 
        in = input('Recompile anyway? (old executable will be copied to the "Old" subdirectory) [n]/y  ','s'); 
        if ~isempty(in) && ( strcmp(in,'y') || strcmp(in,'yes') ) 
            disp('Backing up old executables ...'); 
            if ~exist('Old','file'), mkdir('Old'); end 
            if ~exist(['Old',filesep,'XonOmega.',ext],'file') && exist(filef('XonOmega'),'file') 
                copyfile(filef('XonOmega'),'Old'); 
            end 
            if ~exist(['Old',filesep,'updateSparse.',ext],'file') && exist(filef('updateSparse'),'file') 
                copyfile(filef('updateSparse'),'Old'); 
            end 
            disp('Proceeding to compile...'); 
        else 
            delete(rspfile); 
            return; 
        end 
    end 
     
    disp('Compiling XonOmega.c (file 1 of 2)'); 
    if TESTING 
        if R2006a 
            mex -v -O -DWINDOWS @temp_file.rsp -lmwblas XonOmega.c -output XonOmega 
        else 
            mex -v -O -DWINDOWS @temp_file.rsp -lmwblas XonOmega.c -output XonOmega -largeArrayDims 
        end 
    else 
        if R2006a 
            mex -O -DWINDOWS @temp_file.rsp -lmwblas XonOmega.c -output XonOmega 
        else 
            mex -O -DWINDOWS @temp_file.rsp -lmwblas XonOmega.c -output XonOmega -largeArrayDims 
        end 
    end 
    disp('Success.'); 
    disp('Compiling updateSparse.c (file 2 of 2)'); 
    if TESTING 
        if R2006a 
            mex -v -O -DWINDOWS @temp_file.rsp -lmwblas updateSparse.c -output updateSparse 
        else 
            mex -v -O -DWINDOWS @temp_file.rsp -lmwblas updateSparse.c -output updateSparse -largeArrayDims 
        end 
    else 
        if R2006a 
            mex -O -DWINDOWS @temp_file.rsp -lmwblas updateSparse.c -output updateSparse 
        else 
            mex -O -DWINDOWS @temp_file.rsp -lmwblas updateSparse.c -output updateSparse -largeArrayDims 
        end 
    end 
    disp('Success.'); 
 
    test(); 
     
    disp('Finished installation!'); 
    disp('Note: if you see message about .exp and .lib files not existing, it''s nothing to worry about'); 
     
     
    % and remove rspfile: 
    delete(rspfile); 
 
 
 
% --------------------------------------------------------- 
 
else 
    % tested on 32-bit linux, but hopefully this works 
    % on mac OSX and unix and any 64-bit flavors 
 
 
    % look for necessary files (i.e. make sure we're in right directory) 
    if ~exist(fullfile(pwd,'dbdqr.c'),'file') || ~exist(fullfile(pwd,'bdsqr_mex.c'),'file') 
        disp('Please run this from the same directory as dbdqr.c and bdsqr_mex.c'); 
        return; 
    end 
     
    % look for pre-built executables 
%     mexext('all');  c = computer  
    ext = mexext; 
    filef = @(name) fullfile(pwd,[name,'.',ext]); 
    INSTALL_PROPACK = true; 
    if exist(filef('bdsqr'),'file') 
        disp('Found precompiled PROPACK executables for your platform'); 
        in = input('Recompile anyway? (old executable will be copied to the "Old" subdirectory) [n]/y  ','s'); 
        if ~isempty(in) && ( strcmp(in,'y') || strcmp(in,'yes') ) 
            disp('Backing up old executables ...'); 
            if ~exist('Old','file'), mkdir('Old'); end 
            if ~exist(['Old',filesep,'bdsqr.',ext],'file') 
                copyfile(filef('bdsqr'),'Old'); 
            end 
            disp('Proceeding to compile...'); 
        else 
            in = input('Would you like to install the SVT executables"  [y]/n  ','s'); 
            if ~isempty(in) && ( strcmp(in,'n') || strcmp(in,'no') ) 
                return; 
            end 
            INSTALL_PROPACK = false; 
        end 
    end 
     
    if INSTALL_PROPACK 
        % try to compile, using Matlab's default LAPACK ('libmwlapack') 
        %   but relies on the user having BLAS 
        % Note: we install the fortran versions, since we assume 
        % that a fortran compiler is installed (i.e. included in gcc) 
         
        try 
             
            if NO_FORTRAN 
                % can't compile the reorth.f since this wasn't re-written 
                % in C; dbdqr.f was re-written in c, so compile this. 
                if TESTING 
                    mex -v -O -lmwlapack -lblas bdsqr_mex.c dbdqr.c -DDBDQR_IN_C -output bdsqr 
                else 
                    mex -O -lmwlapack -lblas  bdsqr_mex.c dbdqr.c -DDBDQR_IN_C -output bdsqr 
                end 
            else 
                if TESTING 
                    mex -v -O -lmwlapack -lblas bdsqr_mex.c dbdqr.f -UDBDQR_IN_C -output bdsqr 
                    mex -v -O -lmwlapack -lblas reorth_mex.c reorth.f -output reorth 
                else 
                    mex -O -lmwlapack -lblas  bdsqr_mex.c dbdqr.f -UDBDQR_IN_C -output bdsqr 
                    mex -O -lmwlapack -lblas reorth_mex.c reorth.f -output reorth 
                end 
            end 
         
        catch 
            disp('The problem might be that you don''t have a fortran compiler'); 
            disp('If your system has one but Matlab doesn''t know about it,'); 
            disp('then try running "mex -setup" or directly editing the file'); 
            disp('"mexopt.sh", located in prefdir'); 
            disp('If you still can''t compile with fortran, then you can'); 
            disp('compile bdsqr using the following command:'); 
            disp('  mex -O -lmwlapack -lblas  bdsqr_mex.c dbdqr.c -DDBDQR_IN_C -output bdsqr'); 
            disp('and to deal with reorth, just turn off warnings, as follows:'); 
            disp('  warning(''off'',''PROPACK:NotUsingMex'') '); 
             
            disp('Another way to compile without fortran is to runthe following:'); 
            disp('  global NO_FORTRAN'); 
            disp('  NO_FORTRAN = true;'); 
            disp('and then re-run this script'); 
             
             
            rethrow(lasterror); 
        end 
         
        % PROBLEM: if I link with g77, then symbols have extra _ 
        % So, need to link with gcc or gcc4 
        % Also, order is important: gateway routine should be listed first 
        % Actually, fixing the order means that the g77 v gcc problem is 
        % fixed 
 
        % see if it worked 
        test_PROPACK(); 
    end 
     
    % This is what we don't want to see: PROPACK:NotUsingMex','Using slow 
    % matlab code for bdsqr.' 
     
    % ---------- for SVT ------------ 
    disp('Now compiling SVT mex files'); 
    if exist(filef('XonOmega'),'file') || exist(filef('updateSparse'),'file') 
        disp('Found precompiled SVT executables for your platform'); 
        in = input('Recompile anyway? (old executable will be copied to the "Old" subdirectory) [n]/y  ','s'); 
        if ~isempty(in) && ( strcmp(in,'y') || strcmp(in,'yes') ) 
            disp('Backing up old executables ...'); 
            if ~exist('Old','file'), mkdir('Old'); end 
            if ~exist(['Old',filesep,'XonOmega.',ext],'file') && exist(filef('XonOmega'),'file') 
                copyfile(filef('XonOmega'),'Old'); 
            end 
            if ~exist(['Old',filesep,'updateSparse.',ext],'file') && exist(filef('updateSparse'),'file') 
                copyfile(filef('updateSparse'),'Old'); 
            end 
            disp('Proceeding to compile...'); 
        else 
            return; 
        end 
    end 
     
    disp('Compiling XonOmega.c (file 1 of 2)'); 
    if TESTING 
        if R2006a 
            mex -v -O -lblas XonOmega.c 
        else 
            mex -v -O -lblas XonOmega.c -largeArrayDims 
        end 
    else 
        if R2006a 
            mex -O -lblas XonOmega.c 
        else 
            mex -O -lblas XonOmega.c -largeArrayDims 
        end 
    end 
    disp('Success.'); 
    disp('Compiling updateSparse.c (file 2 of 2)'); 
    if TESTING 
        if R2006a 
            mex -v -O -lblas updateSparse.c 
        else 
            mex -v -O -lblas updateSparse.c -largeArrayDims 
        end 
    else 
        if R2006a 
            mex -O -lblas updateSparse.c 
        else 
            mex -O -lblas updateSparse.c -largeArrayDims 
        end 
    end 
    disp('Success.'); 
 
    test(); 
     
    disp('Finished installation!'); 
     
end 
 
disp('If you would like to install smvp as well, type'); 
disp('  mex -O smvp.c'); 
 
function cc = get_compiler_config() 
    % tested on Windows w/ R2008 only 
    % This has to be in a function, otherwise old versions of matlab 
    % get confused because "mex" is used as a structure (well, a class) 
    % AND as a function. 
    try  
        % this requires both a new verson of matlab and 
        % that a compiler has been selected 
        cc = mex.getCompilerConfigurations('C'); 
        cc = cc.Manufacturer; 
        % Watch out for this error later (in old versions of matlab) 
%         MATLAB:mir_error_function_previously_indexed_by_dot 
    catch 
        cc = []; 
        fprintf('You may want to run ''mex -setup'' to setup the mex compiler,\n if you''ve never used the mex compiler before\n'); 
    end 
 
function test 
 
    disp('Testing XonOmega -- you shouldn''t see any "Warning" messages'); 
    U = randn(10,5); V = randn(10,5); omega = unique( 1+round( 45*rand(15,1) ) ); 
    y1 = XonOmega(U,V,omega); 
    A = U*V'; y2 = A(omega); 
    fprintf('Discrepancy is %f (if < 1e-14, not much to worry about)\n', norm(y1-y2)); 
    disp('Testing XonOmega with complex numbers -- you shouldn''t see any "Warning" messages'); 
    U = randn(10,5) + 1i*randn(10,5); V = randn(10,5) + 1i*randn(10,5); 
    y1 = XonOmega(U,V,omega); A = U*V'; y2 = A(omega); 
    fprintf('Discrepancy is %f (if < 1e-14, not much to worry about)\n', norm(y1-y2)); 
 
    disp('Testing updateSparse -- you shouldn''t see any "Warning" messages'); 
    A = sprand(100,100,.1) + 1i*sprand(100,100,.1);  
    [I,J] = find(A); omega = sub2ind( size(A), I, J ); 
    B = zeros(size(A)); B(omega) = 1;  B = complex(B); 
    B = sparse(B); 
    updateSparse(B,A(omega)); 
    fprintf('Discrepancy is %f (if < 1e-14, not much to worry about)\n', norm(full(B(omega))-A(omega))); 
     
function test_PROPACK 
    disp('Looks like compilation of PROPACK mex file worked.  Now, testing the installation'); 
    a = randn(50,1); b = randn(50,1); 
    [s,bnd]=bdsqr(a,b); 
    A = randn(50); 
    [U,S,V] = svd(A); 
    [U1,S1,V1] = lansvd(A,50,'L'); 
    if norm( diag(S) - diag(S1) ) > 1e-10, disp('Numerical discrepancy!'); end 
    disp('No error messages!  Finished installing PROPACK.');
