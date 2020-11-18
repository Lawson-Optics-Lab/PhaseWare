function PhiUnwrap = UnwrappingConstantin(handles, options)
% PhiUnwrap = cunwrap(Psi, options)
%
% Purpose: Costantini's 2D unwrapping based on Minimum Network Flow
%
% INPUTS:
%   - Psi: 2D-array wrapped phase (radian)
%   - options: optional structure used to modify the behavior of the 
%              unwrapping algorithm. All fields are optional.
%              The fieldname does not need to exactly case-matched.
%       Weight: array of same dimension as the Psi, with is used as
%               the positive weight of the L1 norm of the residual (to be
%               minimized). For example user can provide WEIGHT as function
%               of the interference amplitude or function of phase
%               measurement quality (e.g., fringe contrast).
%               By default: all internal pixels have the same weight of 1,
%               0.5 for edge-pixels, and 0.25 for corner-pixels. 
%       CutSize: even integer scalar, default [4]: the width (in pixel) of
%                the Gaussian kernel use to lower the weight around pixels
%                that do not fulfill rotational relation.
%                If CutSize is 0, no rotational-based weighting will be
%                carried out.
%       RoundK: Boolean, [false] by default. When ROUNDK is true, CUNWRAP
%               forces the partial derivative residuals K1/K2 to be
%               integers.
%       MaxBlockSize: scalar default [125], target linear blocksize.
%                     Set to Inf for single-block global unwrapping.
%                     Note: network flow is costly in CPU for large size.
%       Overlap: scalar default [0.25], fractional overlapping between two
%                neighboring blocks.
%       Verbose: logical [true], control printout information.
%       LPoption: LINPROG option structure as return by OPTIMSET(...).
% OUTPUT:
%   - Phi: 2D-array unwrapped phase (radian)
%
% Minimum network flow computation engine is Matlab LINPROG (optimization
%   toolbox is required)
%
% References:
%   - http://earth.esa.int/workshops/ers97/papers/costantini/
%   - Costantini, M. (1998) A novel phase unwrapping method based on network
%   programming. IEEE Tran. on Geoscience and Remote Sensing, 36, 813-821.
% Author: Bruno Luong <brunoluong@yahoo.com>

Psi = handles;

if nargin<2
    options = struct();
end

if ndims(Psi)>2
    error('CUNWRAP: input Psi must be 2D array.');
end

[ny nx] = size(Psi);

if nx<2 || ny<2
    error('CUNWRAP: size of Psi must be larger than 2');
end

roundK = getoption(options, 'roundk', false);

verbose = getoption(options, 'verbose', true);

% Default weight
w1 = ones(ny,1); w1([1 end]) = 0.5;
w2 = ones(1,nx); w2([1 end]) = 0.5;
weight = w1*w2; % tensorial product
weight = getoption(options, 'weight', weight);

% Split in smaller block size
blocksize = getoption(options, 'maxblocksize', 125);
blocksize = max(blocksize,2);

% Overlapping fraction
p = getoption(options, 'overlap', 0.25);
p = max(min(p,1),0);

% Split the linear index for each dimension
[ix blk_x] = splitidx(blocksize, nx, p);
[iy blk_y] = splitidx(blocksize, ny, p);
nbblk = length(iy)*length(ix);

mydisplay(verbose, 'Number of block(s) = %d\n', nbblk);
mydisplay(verbose, '\tEffective block size = (%d x %d)\n', blk_y, blk_x);
if nbblk>1
    mydisplay(verbose, '\tOverlapping = %2.1f%%\n', p*100);
end

% Allocate arrays
% negative/positive parts of wrapped partial derivatives
x1p = nan(ny-1, nx, class(Psi));
x1m = nan(ny-1, nx, class(Psi));
x2p = nan(ny, nx-1, class(Psi));
x2m = nan(ny, nx-1, class(Psi));

%%
% Loop over the blocks
blknum = 0;
for i=1:length(iy)
    iy0 = iy{i};
    iy1 = iy0(1:end-1);
    iy2 = iy0(1:end);
    
    for j=1:length(ix)
        ix0 = ix{j};
        ix1 = ix0(1:end);
        ix2 = ix0(1:end-1);
        
        blknum = blknum + 1;
        mydisplay(verbose, ['\n' repmat('-', 1, 60) '\n']);
        mydisplay(verbose, 'CUNWRAP: processing block #%d/%d...\n', ...
            blknum, nbblk);
%         set(handles.StatusBox, 'String', ['Processing block ' num2str(blknum) '/' num2str(nbblk)])
        drawnow
        
        options.weight = weight(iy0,ix0);
        options.blknum = blknum;
        % Costantini's minimum network flow resolution for one block
        [x1p(iy1,ix1) x1m(iy1,ix1) x2p(iy2,ix2) x2m(iy2,ix2)] = ...
            cunwrap_blk(Psi(iy0,ix0), ...
            x1p(iy1,ix1), x1m(iy1,ix1),...
            x2p(iy2,ix2), x2m(iy2,ix2),...
            options);
    end
end

%%
% Compute partial derivative Psi1, eqt (1,3)
mydisplay(verbose, '\ty-derivative\n');
i = 1:(ny-1);
j = 1:nx;
[ROW_I ROW_J] = ndgrid(i,j);
I_J = sub2ind(size(Psi),ROW_I,ROW_J);
IP1_J = sub2ind(size(Psi),ROW_I+1,ROW_J);
Psi1 = Psi(IP1_J) - Psi(I_J);
% A priori knowledge that Psi1 is in [-pi,pi)
Psi1 = mod(Psi1+pi,2*pi)-pi;

% Compute partial derivative Psi2, eqt (2,4)
mydisplay(verbose, '\tx-derivative\n');
i = 1:ny;
j = 1:(nx-1);
[ROW_I ROW_J] = ndgrid(i,j);
I_J = sub2ind(size(Psi),ROW_I,ROW_J);
I_JP1 = sub2ind(size(Psi),ROW_I,ROW_J+1);
Psi2 = Psi(I_JP1) - Psi(I_J);
% A priori knowledge that Psi1 is in [-pi,pi)
Psi2 = mod(Psi2+pi,2*pi)-pi;

% Compute the derivative jumps, eqt (20,21)
k1 = x1p-x1m;
k2 = x2p-x2m;

% Round to integer solution (?)
if roundK
    mydisplay(verbose, '\tForce integer 2pi-ambiguities\n');
    k1 = round(k1);
    k2 = round(k2);
end

% Sum the jumps with the wrapped partial derivatives, eqt (10,11)
k1 = reshape(k1,[ny-1 nx]);
k2 = reshape(k2,[ny nx-1]);
k1 = k1+Psi1/(2*pi);
k2 = k2+Psi2/(2*pi);

%%
% Integrate the partial derivatives, eqt (6)
% Not sure why integrate this way is better than otherway around
mydisplay(verbose, 'Integration\n');
% Integration order, not documented
intorder = getoption(options, 'IntOrder', 1);
if intorder==1
    K = cumsum([0 k2(1,:)]);
    K = [K; k1];
    K = cumsum(K,1);
else
    % Integration is this order does not work well (still mysterious for BL)
    K = cumsum([0; k1(:,1)]);
    K = [K k2];
    K = cumsum(K,2);
end

PhiUnwrap = (2*pi)*K;

end % cunwrap

%%
function [x1p x1m x2p x2m] = cunwrap_blk(Psi, ...
                                         X1p, X1m, X2p, X2m, options)
% function [x1p x1m x2p x2m] = cunwrap_blk(Psi, ...
%                                          X1p, X1m, X2p, X2m, options)
% Costantini's minimum network flow resolution for one block

%%
% if getoption(options, 'blknum', NaN) == 8
%     save('debug.mat', 'Psi', 'X1p', 'X1m', 'X2p', 'X2m', 'options');
% end

%%                                    
[ny nx] = size(Psi);

verbose = getoption(options, 'verbose', true);

% the width (in pixel) of the Gaussian kernel to limit effect of
% patch that does not satisfy rotational relation
CutSize = getoption(options, 'cutsize', 4);

% Default weight
w1 = ones(ny,1); w1([1 end])=0.5;
w2 = ones(1,nx); w2([1 end])=0.5;
weight = w1*w2; % tensorial product
weight = getoption(options, 'weight', weight);

% Compute partial derivative Psi1, eqt (1,3)
mydisplay(verbose, '\ty-derivative\n');
i = 1:(ny-1);
j = 1:nx;
[ROW_I ROW_J] = ndgrid(i,j);
I_J = sub2ind(size(Psi),ROW_I,ROW_J);
IP1_J = sub2ind(size(Psi),ROW_I+1,ROW_J);
Psi1 = Psi(IP1_J) - Psi(I_J);
% A priori knowledge that Psi1 is in [-pi,pi)
Psi1 = mod(Psi1+pi,2*pi)-pi;

% Compute partial derivative Psi2, eqt (2,4)
mydisplay(verbose, '\tx-derivative\n');
i = 1:ny;
j = 1:(nx-1);
[ROW_I ROW_J] = ndgrid(i,j);
I_J = sub2ind(size(Psi),ROW_I,ROW_J);
I_JP1 = sub2ind(size(Psi),ROW_I,ROW_J+1);
Psi2 = Psi(I_JP1) - Psi(I_J);
% A priori knowledge that Psi1 is in [-pi,pi)
Psi2 = mod(Psi2+pi,2*pi)-pi;

%%
% The RHS is column-reshaping of a 2D array [ny-1] x [nx-1]
% Build the equality constraint RHS (eqt 17)
mydisplay(verbose, '\tConstraint RHS\n');
beq = 0;
% Compute beq
i = 1:(ny-1);
j = 1:(nx-1);
[ROW_I ROW_J] = ndgrid(i,j);
I_J = sub2ind(size(Psi1),ROW_I,ROW_J);
I_JP1 = sub2ind(size(Psi1),ROW_I,ROW_J+1);
beq = beq + (Psi1(I_JP1)-Psi1(I_J));
I_J = sub2ind(size(Psi2),ROW_I,ROW_J);
IP1_J = sub2ind(size(Psi2),ROW_I+1,ROW_J);
beq = beq - (Psi2(IP1_J)-Psi2(I_J));
beq = -1/(2*pi)*beq;
beq = round(beq);
beq = beq(:);

%%
% The vector of LP is arranged as following:
%   x := (x1p, x1m, x2p, x2m).'
%   x1p, x1m: reshaping of [ny-1] x [nx]
%   x2p, x2m: reshaping of [ny] x [nx-1]
% 
% Row index, used by all foure blocks in Aeq, beq
mydisplay(verbose, '\tConstraint matrix\n');
i = 1:(ny-1);
j = 1:(nx-1);

[ROW_I ROW_J] = ndgrid(i,j);
ROW_I_J = sub2ind([length(i) length(j)],ROW_I,ROW_J);
nS0 = (nx-1)*(ny-1);

% Use by S1p, S1m
[COL_I COL_J] = ndgrid(i,j);
COL_IJ_1 = sub2ind([length(i) length(j)+1],COL_I,COL_J);
[COL_I COL_JP1] = ndgrid(i,j+1);
COL_I_JP1 = sub2ind([length(i) length(j)+1],COL_I,COL_JP1);
nS1 = (nx)*(ny-1);

% Use by S2p, S2m
[COL_I COL_J] = ndgrid(i,j);
COL_IJ_2 = sub2ind([length(i)+1 length(j)],COL_I,COL_J);
[COL_IP1 COL_J] = ndgrid(i+1,j);
COL_IP1_J = sub2ind([length(i)+1 length(j)],COL_IP1,COL_J);
nS2 = (nx-1)*(ny);

% Build four indexes arrays that will be used to split x in four parts
% 29/08/09 bug corrected
offset1p = 0;
idx1p = offset1p+(1:nS1);
offset1m = idx1p(end);
idx1m = offset1m+(1:nS1);
offset2p = idx1m(end);
idx2p = offset2p+(1:nS2);
offset2m = idx2p(end);
idx2m = offset2m+(1:nS2);

% Equality constraint matrix (Aeq)
S1p = + sparse(ROW_I_J, COL_I_JP1,1,nS0,nS1) ...
      - sparse(ROW_I_J, COL_IJ_1,1,nS0,nS1);
S1m = -S1p; 

S2p = - sparse(ROW_I_J, COL_IP1_J,1,nS0,nS2) ...
      + sparse(ROW_I_J, COL_IJ_2,1,nS0,nS2);  
S2m = -S2p;

% Matrix of the LHS of eqt (17)
Aeq = [S1p S1m S2p S2m];
nvars = size(Aeq,2);

% Clean up
clear S1p S1m S2p S2m

%%
% force to be even
CutSize = ceil(CutSize/2)*2; 

% Modify weight to limit the effect of points that violate
% the rorational condition. The weight is graduataly increase
% around the points.
if CutSize>0
    mydisplay(verbose, '\tAdjust weight (CutSize = %d pixels)\n', CutSize);
    % Truncated Gaussian Kernel
    v = 1*linspace(-1,1,CutSize);
    [x y] = meshgrid(v,v);
    kernel = 1.1*exp(-(x.^2+y.^2));

    rotdegradation = double(reshape(beq~=0, [ny nx]-1)); % 0 or 1
    mydisplay(verbose, '\tInconsistency = %d/%d\n', ...
              sum(rotdegradation(:)), numel(beq));
    rotdegradation = conv2(rotdegradation, kernel, 'full');
    rotdegradation = rotdegradation(CutSize/2 + (0:ny-1),...
                                    CutSize/2 + (0:nx-1));
    % mininum weight threshold      
    wmin = 1e-2; % weight never goes under this value, small positive value
    rotdegradation = min(rotdegradation,1-wmin);
    weight = weight .* (1-rotdegradation);
end

c1 = 0.5*(weight(1:ny-1,:)+weight(1:ny-1,:));
c2 = 0.5*(weight(:,1:nx-1)+weight(:,1:nx-1));

% Cost vector, eqt (16)
cost = zeros(nvars,1);
cost(idx1p) = c1(:);
cost(idx1m) = c1(:);
cost(idx2p) = c2(:);
cost(idx2m) = c2(:);

%%
% Lower and upper bound, eqt (18,19)
L = zeros(nvars,1);
U = Inf(size(L)); % No upper bound, U=[];

%%
% Lock the x values to prior computed value (from calculation on other
% blocks)
i1 = find(~isnan(X1p));
i2 = find(~isnan(X2p));
mydisplay(verbose, '\tLock variables = %d/%d\n', ...
          2*(length(i1)+length(i2)), size(Aeq,2));
      
% Lock method, not documented
lockadd = 1;
lockremove = 2; %#ok
lockmethod = getoption(options, 'lockmethod', lockadd);

if lockmethod==lockadd
    % Lock matrix and values, eqt (26, 27), larger system
    L1p = sparse(1:length(i1), idx1p(i1), 1, length(i1), size(Aeq,2));
    L1m = sparse(1:length(i1), idx1m(i1), 1, length(i1), size(Aeq,2));
    L2p = sparse(1:length(i2), idx2p(i2), 1, length(i2), size(Aeq,2));
    L2m = sparse(1:length(i2), idx2m(i2), 1, length(i2), size(Aeq,2));
    
    AL = [L1p; L1m; L2p; L2m];
    bL = [X1p(i1); X1m(i1); X2p(i2); X2m(i2)];
    clear L1p L1m L2p L2m % clean up
    
    % Find the rows in Aeq with all variables locked
    ColPatch = [offset1p+COL_IJ_1(:) offset1p+COL_I_JP1(:) ...
                offset2p+COL_IJ_2(:) offset2p+COL_IP1_J(:)];
    
    [trash ColDone] = find(AL); %#ok
    remove = all(ismember(ColPatch, ColDone),2);
    Aeq = [Aeq(~remove,:); AL];
    beq = [beq(~remove,:); bL];
    
    % No need to bother with what already computed
    cost(idx1p(i1)) = 0;
    cost(idx1m(i1)) = 0;
    cost(idx2p(i2)) = 0;
    cost(idx2m(i2)) = 0;
    
    L(idx1p(i1)) = -Inf;
    L(idx1m(i1)) = -Inf;
    L(idx2p(i2)) = -Inf;
    L(idx2m(i2)) = -Inf;
    
    clear AL bL trash ColPatch ColDone remove i1 i2 % clean up
else
    % Lock by remove the overlapped variables, smaller system
    % But *seems* more affected by error propagation
    % BL think both method should be strictly equivalent (!)
    
    % remove the equality contribution from the RHS
    lock = zeros(nvars,1,class(Aeq));
    lock(idx1p(i1)) = X1p(i1);
    lock(idx1m(i1)) = X1m(i1);
    lock(idx2p(i2)) = X2p(i2);
    lock(idx2m(i2)) = X2m(i2);
    beq = beq - Aeq*lock;
    
    % Remove the variables
    vremove = [idx1p(i1) idx1m(i1) idx2p(i2) idx2m(i2)];
    % keep is use later to dispatch partial derivative
    keep = true(nvars,1); keep(vremove) = false;
    Aeq(:,vremove) = [];
    L(vremove) = []; U(vremove) = [];
    cost(vremove) = [];
    
    % Find the rows in Aeq with all variables locked
    ColPatch = [offset1p+COL_IJ_1(:) offset1p+COL_I_JP1(:) ...
                offset2p+COL_IJ_2(:) offset2p+COL_IP1_J(:)];
    
    % Remove the equations
    eremove = all(ismember(ColPatch, vremove),2);
    Aeq(eremove,:) = [];
    beq(eremove,:) = [];
    
    clear vremove eremove ColPatch lock % clean up
end

%%
% To do: implement Bertsekas/Tseng's relaxation method, ref. [20]
% Call LP solver
mydisplay(verbose, '\tMinimum network flow resolution\n');
mydisplay(verbose, '\t\tmatrix size = (%d,%d)\n', size(Aeq,1),size(Aeq,2));
if ~isempty(which('linprog'))
    % Call Matlab Linprog
    % Adjust optional LP options at your preference, see
    % http://www.mathworks.com/access/helpdesk/help/toolbox/optim/ug/linprog.html
    % http://www.mathworks.com/access/helpdesk/help/toolbox/optim/ug/optimset.html
    
    mydisplay(verbose, '\tLINPROG...\n');
    
    % BL has not checked because he does not have the right toolbox
    LPoption = getoption(options, 'LPoption', {});
    if ~iscell(LPoption)
        LPoption = {LPoption};
    end
    % LPoption = {optimset(...)}
    sol = linprog(cost,[],[],Aeq,beq,L,U,[],LPoption{:});
elseif ~isempty(which('BuildMPS'))
    
    % Here is BL Linprog, call Matlab linprog instead to get "SOL",
    % the solution of the above LP problem
    mpsfile='costantini.mps';
    
    mydisplay(verbose, '\tConversion to MPS file\n');
    Contain = BuildMPS([], [], Aeq, beq, cost, L, U, 'costantini');
    OK = SaveMPS(mpsfile,Contain);
    if ~OK
        error('CUNWRAP: Cannot write mps file');
    end
    
    PCxpath=App_path('PCx');
    [OK outfile]=invokePCx(mpsfile,PCxpath,verbose==0);
    if ~OK
        mydisplay(verbose, 'PCx path=%s\n', PCxpath);
        mydisplay(verbose, 'Cannot invoke LP solver, PCx not installed?\n');
        error('CUNWRAP: Cannot invoke PCx');
    end
    [OK sol]=readPCxoutput(outfile);
    if ~OK
        error('CUNWRAP: Cannot read PCx outfile, L1 fit might fails.');
    end
else
    error('CUNWRAP: cannot detect network flow (LP) engine');
end

%%
% Displatch the LP solution
if lockmethod==lockadd
    x = sol;
else
    x = zeros(size(keep),class(sol));
    x(keep) = sol;
end
x1p = reshape(x(idx1p), [ny-1 nx]);
x1m = reshape(x(idx1m), [ny-1 nx]);
x2p = reshape(x(idx2p), [ny nx-1]);
x2m = reshape(x(idx2m), [ny nx-1]);

end

%% Get defaut option
function value = getoption(options, name, defaultvalue)
% function value = getoption(options, name, defaultvalue)
    value = defaultvalue;
    fields = fieldnames(options);
    found = strcmpi(name,fields);
    if any(found)
        i = find(found,1,'first');
        if ~isempty(options.(fields{i}))
            value = options.(fields{i});
        end
    end
end

%%
function [ilist blocksize] = splitidx(blocksize, n, p)
% function ilist = splitidx(blocksize, n, p)
%
% return the cell array, each element is subindex of (1:n)
% The union is (1:n) with overlapping fraction is p (0<p<1)

if blocksize>=n
    ilist = {1:n};
    blocksize = n;
else
    q = 1-p;

    % Number of blocks
    k = (n/blocksize - p) / q;
    k = ceil(k);

    % Readjust the block size, float
    blocksize = n/(k*q + p);

    % first index
    firstidx = round(linspace(1,n-ceil(blocksize)+1, k));
    lastidx = round(firstidx+blocksize-1);
    lastidx(end) = n;
    % Make sure they are overlapped
    lastidx(1:end-1) = max(lastidx(1:end-1),firstidx(2:end));

    % Put the indexes of k blocks into cell array
    ilist = cell(1,length(firstidx));
    for k=1:length(ilist)
        ilist{k} = firstidx(k):lastidx(k);
    end
    
    blocksize = round(blocksize);
end
end

% function to mydisplay a text on command line window
function mydisplay(verbose, varargin)
if verbose
    fprintf(varargin{:});
end
end

