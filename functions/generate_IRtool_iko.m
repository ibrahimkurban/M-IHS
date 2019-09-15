function [A, b0, x0, probinfo] = generate_IRtool_iko(N, problem, example, varargin)
%IRtool generates Hansen data by given progile number and varargin
% parameters
%
%    [A, b0, x0, probinfo] = generate_IRtool(n, problem, example, varargin)
%
% VARARGIN
% 'type', double
% 'level', double
%
% PROBLEMS
%   1 : BLURRING  
%           example : 1,...,9
%           'type', : 1,...,6
%           'level' : 1,2,3 (blurring level)
%   2 : SEISMIC  
%           example : 1,...,8
%           'type', : 1,2
%           'level' : 1,2,3... makes overdetermined
%   3 : X-RAY TOMOGRAPHIC  
%           example : 1,...,8
%           'type', : 1,2
%           no level input
%   4 : SPHERICAL  
%           example : 1,...,8
%           no type  input
%           no level input
%
%

switch problem
    case 1
        %%
        fprintf('BLURRING PROBLEM is chosen\n')
        examples    = {'pattern1', 'pattern2', 'sppattern', 'ppower', 'smooth', 'dot2', 'dotk', 'satellite', 'hst'};
        types       = {'gauss', 'defocus', 'speckle', 'shake', 'motion', 'rotation'};
        level       = {'mild', 'medium', 'severe'};
        
        if( nargin-length(varargin)+1 < 3 )
            exa = 'pattern1';
            psf = 'gauss';
            lev = 'mild';
        elseif(nargin-length(varargin)+1 == 3)
            exa = examples{example};
            psf = 'gauss';
            lev = 'mild';
        else
            exa = examples{example};
            psf = 'gauss';
            lev = 'mild';
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case 'type'
                        psf = types{varargin{i+1}};
                    case 'level'
                        lev = level{varargin{i+1}};
                    otherwise
                        fprintf('\n!!!!!!!!!! ILLEGAL INPUT !!!!!!!!!!\n\n')
                end
            end
        end
        %options
        options = PRset('trueImage', exa, 'PSF', psf, 'BlurLevel', lev, ...
            'CommitCrime', 'on');
        [A, b0, x0, probinfo] = PRblur(N, options);
        
        if(~issparse(A))
            fprintf('\n!!!!!!!!!! WARNING: DATA IS NOT SPARSE !!!!!!!!!!\n\n')
            A      = sparse(full(A));
        end
        [n,d] = size(A);
        fprintf('---- Example       : %s\n', exa);
        fprintf('---- PSF           : %s\n', psf);
        fprintf('---- Blur Level    : %s\n', lev);
        fprintf('---- Commit Crime  : on\n');
        fprintf('---- Boundary      : reflective\n');
        fprintf('---- Dimensions    : %1.2e x %1.2e\n\n', n,d);
        
    case 2
        %%
        fprintf('SEISMIC TRAVEL-TIME TOMOGRAPHY PROBLEM is chosen\n')
        examples    = {'tectonic', 'smooth', 'binary', 'threephases', 'threephasessmooth', 'fourphases', 'grains', 'ppower'};
        types       = {'ray', 'fresnel'};
        addit       = {};
        if( nargin-length(varargin)+1 < 3 )
            exa = 'tectonic';
            psf = 'ray';
            lev = 1;
        elseif(nargin-length(varargin)+1 == 3)
            exa = examples{example};
            psf = 'ray';
            lev = 1;
        else
            exa = examples{example};
            psf = 'ray';
            lev = 1;
            k       = -1;
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case 'type'
                        psf = types{varargin{i+1}};
                    case 'level'
                        lev = varargin{i+1};
                    otherwise
                        k            = k+2;
                        addit(k:k+1) = varargin(i:i+1); 
                end
            end
        end
        %options
        options = PRset('phantomImage', exa, 'wavemodel', psf, 's', lev*N, addit{:});
        [A, b0, x0, probinfo] = PRseismic(N, options);
        
        [n,d] = size(A);
        fprintf('---- Example       : %s\n', exa);
        fprintf('---- Wavemodel     : %s\n', psf);
        fprintf('---- # of sources  : %d*n\n', lev);
        fprintf('-----# of receiver : %d*n\n', options.p/N)
        fprintf('---- Dimensions    : %1.2e x %1.2e\n\n', n,d);
        
    case 3
        %%
        fprintf('X-RAY TOMOGRAPHIC PROBLEM is chosen\n')
        examples    = {'shepplogan', 'smooth', 'binary', 'threephases', 'threephasessmooth', 'fourphases', 'grains', 'ppower'};
        types       = {'parallel', 'fancurved'};
        addit       = {};
        if( nargin-length(varargin)+1 < 3 )
            exa = 'shepplogan';
            psf = 'parallel';
            
        elseif(nargin-length(varargin)+1 == 3)
            exa = examples{example};
            psf = 'parallel';
        else
            k   = -1;
            exa = examples{example};
            psf = 'parallel';
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case 'type'
                        psf = types{varargin{i+1}};
                    otherwise
                        k            = k+2;
                        addit(k:k+1) = varargin(i:i+1); 
                end
            end
            
        end
        %options
        options = PRset('phantomImage', exa, 'CTtype', psf, addit{:});
        [A, b0, x0, probinfo] = PRtomo(N, options);
        
        [n,d] = size(A);
        fprintf('---- Example       : %s\n', exa);
        fprintf('---- CTtype        : %s\n', psf);
        fprintf('---- Dimensions    : %1.2e x %1.2e\n\n', n,d);
        
    case 4
        %%
        fprintf('SPHERICAL MEANS TOMOGRAPHY PROBLEM is chosen\n')
        examples    = {'shepplogan', 'smooth', 'binary', 'threephases', 'threephasessmooth', 'fourphases', 'grains', 'ppower'};
        
        if( nargin-length(varargin) < 3 )
            exa = 'shepplogan';
        elseif(nargin-length(varargin) == 3)
            exa = examples{example};
        else
            exa = examples{example};
            fprintf('\n!!!!!!!!!! ILLEGAL INPUT !!!!!!!!!!\n\n')
        end
        %options
        options = PRset('phantomImage', exa);
        [A, b0, x0, probinfo] = PRspherical(N, options);
        
        [n,d] = size(A);
        fprintf('---- Example       : %s\n', exa);
        fprintf('---- Dimensions    : %1.2e x %1.2e\n\n', n,d);
        
end
b0  = A*x0;
end





