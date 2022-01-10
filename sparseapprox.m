function varargout = sparseapprox(X, D, met, varargin)

% The coefficients or weights, W, are usually (but not always) sparse, 
% i.e. number of non-zero coefficients below a limit,
% some methods use the 1-norm but these are not much tested.
% 
% The reconstrucion or approximation of X is (D*W).
% The approximation error is R = X - D*W;

mfile = 'sparseapprox';

if size(X,2)~=64
    I = randperm(size(X,2));
    X_test_sparsity=X;
    X=(X(:,I(1:64))); 
end

%% defaults, initial values
tstart = tic;

    
[N,L] = size(X);
K = size(D,2);
norm2X = sqrt(sum(X.*X)); % ||x(i)||_2     1xL
W = zeros(K,L);      % the weights (coefficients)
tnz = ceil(N/2)*ones(1,L); % target number of non-zeros
thrActive = false;    % is set to true if tnz, tre or tae is given
                      % and used for methods: pinv, backslash, linprog and
                      % FOCUSS
doOP = true;          % do Orthogonal Projection when thresholding                 
relLim = 1e-6;
tre = relLim*ones(1,L);   % target relative error: ||r|| <= tre*||x||

verbose = 0;
done = false;
    
%% get the options
nofOptions = nargin-3;
optionNumber = 1;
fieldNumber = 1;
while (optionNumber <= nofOptions)
    if isstruct(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        sNames = fieldnames(sOptions);
        opName = sNames{fieldNumber};
        opVal = sOptions.(opName);
        % next option is next field or next (pair of) arguments
        fieldNumber = fieldNumber + 1;  % next field
        if (fieldNumber > numel(sNames)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    elseif iscell(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        opName = sOptions{fieldNumber};
        opVal = sOptions{fieldNumber+1};
        % next option is next pair in cell or next (pair of) arguments
        fieldNumber = fieldNumber + 2;  % next pair in cell
        if (fieldNumber > numel(sOptions)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    else  
            opName = varargin{optionNumber};
            opVal = varargin{optionNumber+1};
            optionNumber = optionNumber + 2;  % next pair of options
        
    end
    % interpret opName and opVal
    if strcmpi(opName,'targetNonZeros') || strcmpi(opName,'tnz')
        if strcmpi(met,'GMP') 
            tnz = opVal;   % GMP will distribute the non-zeros
        else
            if numel(opVal)==1
                tnz = opVal*ones(1,L);
            elseif numel(opVal)==L
                tnz = reshape(opVal,1,L);
            else
                error([mfile,': illegal size of value for option ',opName]);
            end
        end
        thrActive = true;
    end
    if strcmpi(opName,'targetRelativeError') || strcmpi(opName,'tre')
        if numel(opVal)==1
           tre = opVal*ones(1,L);
        elseif numel(opVal)==L
           tre = reshape(opVal,1,L);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
        thrActive = true;
    end
    if strcmpi(opName,'targetAbsoluteError') || strcmpi(opName,'tae')
        if numel(opVal)==1
           tae = opVal*ones(1,L);
        elseif numel(opVal)==L
           tae = reshape(opVal,1,L);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
        thrActive = true;
    end
    if ( strcmpi(opName,'nIt') || strcmpi(opName,'nofIt') || ...
         strcmpi(opName,'numberOfIterations') )
        if (isnumeric(opVal) && numel(opVal)==1)
            nIt = max(floor(opVal), 1);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'p') || strcmpi(opName,'pFOCUSS')
        if (isnumeric(opVal) && numel(opVal)==1)
            pFOCUSS = min(opVal, 1);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'l') || strcmpi(opName,'lambdaFOCUSS')
        if (isnumeric(opVal) && numel(opVal)==1)
            lambdaFOCUSS = abs(opVal);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'nComb')
        if (isnumeric(opVal) && numel(opVal)==1)
            nComb = max(floor(opVal), 2);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'paramSPAMS')
        if (isstruct(opVal))
            paramSPAMS = opVal;
        else
            error([mfile,': option paramSPAMS is not a struct as it should be, see SPAMS help.']);
        end
    end
    if strcmpi(opName,'globalReDist')
        if (isnumeric(opVal) && numel(opVal)==1)
            globalReDist = min(max(floor(opVal), 0), 2);  % 0, 1 or 2
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'doOP') 
        if (islogical(opVal)); doOP = opVal; end;
        if isnumeric(opVal); doOP = (opVal ~= 0); end;
    end
    if strcmpi(opName,'GMPLoopSize')
        if (isnumeric(opVal) && numel(opVal)==1)
            GMPLoopSize = max(floor(opVal), 2);
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'tSSE') || strcmpi(opName,'targetSSE')
        if (isnumeric(opVal) && numel(opVal)==1)
            targetSSE = min(max(opVal, 0), sum(sum(X.*X)));
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'tSNR') || strcmpi(opName,'targetSNR')
        if (isnumeric(opVal) && numel(opVal)==1)
            targetSSE = 10^(-abs(opVal)/10) * sum(sum(X.*X));
        else
            error([mfile,': illegal size of value for option ',opName]);
        end
    end
    if strcmpi(opName,'verbose') || strcmpi(opName,'v')
        if (islogical(opVal) && opVal); verbose = 1; end;
        if isnumeric(opVal); verbose = opVal(1); end;
    end
end

if exist('tae','var')   % if both exist 'tae' overrules 'tre'
    tre = tae./norm2X;
elseif exist('tre','var')  % 'tre' was given a default value
    tae = tre.*norm2X;
else                       % so this case is redundant
    disp(' ??? This is never printed.');
end

%% Display info
if (verbose > 1)  % very verbose
    disp(' ');
    disp([mfile,' with method ',met,' started']);
    
end

%%  Method of Dictionary
if strcmpi(met,'MOF') || strcmpi(met,'MethodOfFrames') || strcmpi(met,'Graph Regulized')
    textMethod = 'Method of Graph Regulized';
    if (verbose >= 1); disp([mfile,': ',textMethod]); end;
    W = pinv(D)*X;
    if thrActive  % then adjust w by setting more to zero
        W = setSmallWeightsToZero(X,D,W,tnz,tae,doOP);
    end
    done = true;
end

if  strcmpi(met,'Test')
     W=X_test_sparsity(size(X_test_sparsity,1),size(X_test_sparsity,2),1);
end
%% Now we are finished, 'done' should be true
%  but test this before finding (and/or) displaying results
% W = full(W);   

%% may display info before returning
if done && ((verbose > 1) || (nargout >= 2))  % need some results
    R = X - D*W;
    varX = var(X(:));     
    varR = var(R(:));     
    if (varR > 0)
        snr = 10*log10(varX/varR);
    else
        snr = inf;
    end
    norm0X = sum(X ~= 0);
    norm1X = sum(abs(X));
    normiX = max(abs(X));  
    norm0R = sum(R ~= 0);
    norm1R = sum(abs(R));
    norm2R = sqrt(sum(R.*R));  
    normiR = max(abs(R));  
    norm0W = sum(W ~= 0);
    norm1W = sum(abs(W));
    norm2W = sqrt(sum(W.*W));  
    normiW = max(abs(W));  
end
if done && (verbose >= 2)  % very verbose
    if (snr < 100)
        disp([mfile,': SNR for the reconstruction is ',...
              num2str(snr,'%7.4f')]);
    elseif (snr < 500)
        disp([mfile,': Almost perfect reconstruction, SNR > 100.']);
    else
        disp([mfile,': Perfect reconstruction, X = D*W.']);
    end
    disp(['Number of non-zeros in W is ',int2str(sum(norm0W)),...
        ' (sparseness factor is ',num2str(sum(norm0W)/(N*L)),')']);
    if exist('numberOfIterations','var');        
        disp(['Average number of iterations for each column : ',...
            num2str(mean(numberOfIterations),'%5.1f')]);
    end
    %
    disp(['X: ',num2str(min(norm0X)),' <= ||x||_0 <= ',...
        num2str(max(norm0X)),'  and mean is ',num2str(mean(norm0X))]);
    disp(['   ',num2str(min(norm1X)),' <= ||x||_1 <= ',...
        num2str(max(norm1X)),'  and mean is ',num2str(mean(norm1X))]);
    disp(['   ',num2str(min(norm2X)),' <= ||x||_2 <= ',...
        num2str(max(norm2X)),'  and mean is ',num2str(mean(norm2X))]);
    disp(['   ',num2str(min(normiX)),' <= ||x||_inf <= ',...
        num2str(max(normiX)),'  and mean is ',num2str(mean(normiX))]);
    disp(['R: ',num2str(min(norm0R)),' <= ||r||_0 <= ',...
        num2str(max(norm0R)),'  and mean is ',num2str(mean(norm0R))]);
    disp(['   ',num2str(min(norm1R)),' <= ||r||_1 <= ',...
        num2str(max(norm1R)),'  and mean is ',num2str(mean(norm1R))]);
    disp(['   ',num2str(min(norm2R)),' <= ||r||_2 <= ',...
        num2str(max(norm2R)),'  and mean is ',num2str(mean(norm2R))]);
    disp(['   ',num2str(min(normiR)),' <= ||r||_inf <= ',...
        num2str(max(normiR)),'  and mean is ',num2str(mean(normiR))]);
    disp(['W: ',num2str(min(norm0W)),' <= ||w||_0 <= ',...
        num2str(max(norm0W)),'  and mean is ',num2str(mean(norm0W))]);
    disp(['   ',num2str(min(norm1W)),' <= ||w||_1 <= ',...
        num2str(max(norm1W)),'  and mean is ',num2str(mean(norm1W))]);
    disp(['   ',num2str(min(norm2W)),' <= ||w||_2 <= ',...
        num2str(max(norm2W)),'  and mean is ',num2str(mean(norm2W))]);
    disp(['   ',num2str(min(normiW)),' <= ||w||_inf <= ',...
        num2str(max(normiW)),'  and mean is ',num2str(mean(normiW))]);
    disp(' ');
    disp([mfile,' with ',met,' done. Used time is ',num2str(toc(tstart))]);
    disp(' ');
end

%% Return Outputs
if done
    if (nargout >= 1)
        varargout{1} = W;
    end
    if (nargout >= 2)
        varargout{2} = struct( 'time', toc(tstart), 'snr', snr, ...
            'textMethod', textMethod, ...
            'norm0X', norm0X, 'norm1X', norm1X, ...
            'norm2X', norm2X, 'normiX', normiX, ...
            'norm0R', norm0R, 'norm1R', norm1R, ...
            'norm2R', norm2R, 'normiR', normiR, ...
            'norm0W', norm0W, 'norm1W', norm1W, ...
            'norm2W', norm2W, 'normiW', normiW );
        % extra output arguments for standard and regularized FOUCUSS 
        if exist('sparseInW','var');        
            varargout{2}.sparseInW = sparseInW;
            varargout{2}.edges = edges;
        end
        if exist('changeInW','var');        
            varargout{2}.changeInW = changeInW;
        end
        if exist('numberOfIterations','var');        
            varargout{2}.numberOfIterations = numberOfIterations;
        end
        if exist('sol','var');        
            varargout{2}.sol = sol;
        end
        if exist('paramSPAMS','var');        
            varargout{2}.paramSPAMS = paramSPAMS;
        end
    end
else  %  ~done
    
    if  strcmpi(met,'Test')
        switch W
            case 1
                W='sibe_salem_ghermez';
            case 2
                W='sibe_salem_zard';            
            case 3
                W='gal_va_kerme_sibe_zard';
            case 4
                W='gal_va_senake_sibe_ghermez';
            case 5
                W='gal_va_senake_sibe_zard';
            case 6
                W='gal_va_zangare_sibe_ghermez';
            case 7
                W='gal_va_zangare_sibe_zard';
            case 8
                W='gale_ghermez';
            case 9
                W='gale_zard';
            case 10
                W='kerm_va_senake_sibe_ghermez';
            case 11
                W='kerm_va_senake_sibe_zard';
            case 12
                W='kerm_va_zangare_sibe_zard';
            case 13
                W='kerme_ghermez';
            case 14
                W='kerme_zard';            
            case 15
                W='senak_va_zangare_sibe_zard';
            case 16
                W='senake_ghermez';
            case 17
                W='senake_zard';
            case 18
                W='zangare_zard';
            otherwise
                disp('Image is not preprocessed, Please do the filters on it');
                W='???';
        end
        varargout{1} = W;         
    else
        textMethod = [met,'. Coefficients W not found.'];
        disp([mfile,': ',textMethod]);
        if (nargout >= 1)
            varargout{1} = [];
        end
        if (nargout >= 2)
            varargout{2} = struct( 'time', toc(tstart), ...
                'textMethod', textMethod );
        end
    end
end
return

%%
function W = setSmallWeightsToZero(X,D,W,tnz,tae,doOP)
[K,L] = size(W);  % K=size(D,2), L=size(X,2)
for i = 1:L
    [absw,I] = sort(abs(W(:,i)),'descend');
    w = zeros(K,1);
    for s = 1:tnz(i)
        if (absw(s) < 1e-10); break; end;
        if doOP % do orthogonal projection onto the selected columns of D
            Di = D(:,I(1:s));
            w(I(1:s)) = (Di'*Di)\(Di'*X(:,i));
        else  % simple thresholding keeps the largest coefficients unaltered
            w(I(s)) = W(I(s),i);
        end
        r = X(:,i) - D*w;
        if (sqrt(r'*r) < tae(i)); break; end;
    end
    W(:,i) = w;
end
return
        

