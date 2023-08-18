function NCA_estimate_stat(job)
% Main NCA function for estimating NCA statistics
%
% Existing SPM designs (saved as SPM.mat) are used to recognise the model 
% and restrict the NCA to regression designs. Additional nuisance parameters are
% ignored, as the necessity of a factor does not depend on other factors, 
% and the model specification can be restricted to the potentially necessary 
% conditions of interest. Therefore, we do not need to control for other 
% variables because the estimated effect size of the necessity variable
% is not contaminated by the presence or absence of other (necessary) 
% variables. 
%
% FORMAT NCA_estimate_stat(job)
% job        - job from tbx_cfg_NCA
%
% job fields:
%   data
%
%   nproc
%     0 - no multi-threading
%     x - number of processors for multi-threading
%
%   conspec.contrasts
%     Inf - interactive selection of contrast
%     x   - index of contrast(s)
%
%   mask
%     mask for restricting NCA estimation (SVC)
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% disable parallel processing for only one SPM.mat file
if numel(job.data) == 1
  job.nproc = 0;
end

% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
  cat_parallelize(job,mfilename,'data');
  return
% or run though all data step-by-step
elseif numel(job.data)>1
  data = job.data;
  for i=1:numel(job.data)
    job = rmfield(job,'data');
    job.process_index = i;
    job.data{1} = data{i};
    NCA_estimate_stat(job)
  end
  return
end

% Use tail approximation from the Gamma distribution for corrected P-values 
use_gamma_tail_approximation = true;
tail_approximation_wo_unpermuted_data = false;

% display permuted design matrix (otherwise show t distribution)
show_permuted_designmatrix = true;
    
% colors and alpha levels
col   = [0.25 0 0; 1 0 0; 1 0.75 0];
alpha = [0.05      0.01   0.001];
    
% give same results each time
if exist('rng','file') == 2
  rng('default')
  rng(0)
else
  rand('state',0);
end

% tolerance for comparing real numbers
tol = 1e-4;
    
load(job.data{1});
cwd = fileparts(job.data{1});

%-Check that model has been estimated
if ~isfield(SPM, 'xVol')
  str = { 'This model has not been estimated.';...
              'Would you like to estimate it now?'};
  if spm_input(str,1,'bd','yes|no',[1,0],1)
    cd(cwd)
    SPM = spm_spm(SPM);
  else
    return
  end
end
    
Ic0 = job.conspec.contrasts;

% if just one contrast is defined use this and skip interactive selection
if (numel(Ic0) == 1) && ~isfinite(Ic0) && isfield(SPM,'xCon') && (numel(SPM.xCon) == 1)
  Ic0 = 1;
end

% check whether contrast are defined
if (numel(Ic0) == 1)
  if ~isfinite(Ic0) || ~isfield(SPM,'xCon') || (isfield(SPM,'xCon') && isempty(SPM.xCon))
    [Ic0,xCon] = spm_conman(SPM,'T',Inf,...
          '  Select contrast(s)...',' ',1);
    SPM.xCon = xCon;
  end
end

% for default SPM.mat results has to be called first
if isempty(SPM.xCon(Ic0(1)).eidf)
  fprintf('You have to call results first.\n');
  cat_spm_results_ui('Setup',SPM);
  load(job.data{1});
end

% check that no temporal filter was used
if isstruct(SPM.xX.K)
  fprintf('ERROR: No first level analysis with temporal correlations allowed.\n');
  return
end

% get some parameters from SPM
xX     = SPM.xX;
VY     = SPM.xY.VY;
n_data = size(xX.X,1);

% sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design
xX = correct_xX(xX);

% check for longitudinal designs (paired t-test, flexible factorial)
repeated_anova = ~isempty(xX.iB);
if repeated_anova
  error('NCA cannot be applied for repeated measure Anova.')
else
  exch_block_labels = ones(1,n_data);
end

% check for meshes
if spm_mesh_detect(VY)
  mesh_detected = true;
else
  mesh_detected = false;
end

% check for mask image that should exist for any analysis
if exist(fullfile(cwd, 'mask.img'),'file')
  file_ext = '.img';
elseif exist(fullfile(cwd, 'mask.nii'),'file')
  file_ext = '.nii';
elseif exist(fullfile(cwd, 'mask.gii'),'file')
  file_ext = '.gii';
else
  fprintf('WARNING: No mask file found. Switch to test mode.\n');
end

% check whether CAT12 is installed
if mesh_detected
  if ~exist('spm_cat12','file')
    error('For using surface analysis you need to install CAT12.');
  end
  
end

% get mask file
if isempty(job.mask)
  maskname = fullfile(cwd,['mask' file_ext]);
else
  if ~isempty(job.mask{1})
    maskname = job.mask{1};
  else    
    maskname = fullfile(cwd,['mask' file_ext]);
  end
end

% load mask
try
  Vmask = spm_data_hdr_read(maskname);
catch
  if mesh_detected
    maskname = spm_select(1,'mesh','select surface mask');
  else
    maskname = spm_select(1,'image','select mask image');
  end
  Vmask = spm_data_hdr_read(maskname);
end
  
% if first image was not found you have to select all files again
if ~exist(VY(1).fname,'file')

  fprintf('Data not found. Please select data in the order defined in the design matrix.\n');
  n = size(SPM.xY.VY,1);
  if mesh_detected
    P = spm_select(n,'mesh','select surfaces');
  else
    P = spm_select(n,'image','select images');
  end
  
  VY = spm_data_hdr_read(P);
  
  %-Apply gSF to memory-mapped scalefactors to implement scaling
  %--------------------------------------------------------------------------
  for i = 1:n
    VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*SPM.xGX.gSF(i);
    if mesh_detected
      if isfield(VY(i).private.private.data{1}.data,'scl_slope')
        VY(i).private.private.data{1}.data.scl_slope = ...
            VY(i).private.private.data{1}.data.scl_slope * SPM.xGX.gSF(i);
        VY(i).private.private.data{1}.data.scl_inter = ...
            VY(i).private.private.data{1}.data.scl_inter * SPM.xGX.gSF(i);
      end
    else
      VY(i).private.dat.scl_slope = ...
          VY(i).private.dat.scl_slope * SPM.xGX.gSF(i);
      VY(i).private.dat.scl_inter = ...
          VY(i).private.dat.scl_inter * SPM.xGX.gSF(i);
    end
  end
  
  SPM.xY.VY = VY;
    
  % update SPM
  if size(SPM.xY.VY,1)==n
    save(job.data{1},'SPM','-v7.3');
  else
    fprintf('ERROR: Number of files is not correct\n');
    return
  end
end
      
% check whether mask images fits to the data
if mesh_detected, dim_index = 1; else dim_index=1:3; end
if sum(sum((Vmask.mat-VY(1).mat).^2)) > 1e-6 || any(Vmask.dim(dim_index) ~= VY(1).dim(dim_index))
  error('Mask must have the same dimensions and orientation as the data.');
end

% read mask and data
mask = spm_data_read(Vmask);

ind_mask = find(mask>0);
n = numel(VY);

if ~isempty(ind_mask)
  Y = zeros([length(ind_mask) n],'single');
  C = [];
  
  % load data
  fprintf('Load data ')
  for i=1:n
    fprintf('.')
    tmp = spm_data_read(VY(i));
    Y(:,i) = tmp(ind_mask);
  end
  
  fprintf('\n')
  clear tmp;

  % whitening matrix
  W = single(full(xX.W));
  Y = Y*W;

else
  error('Empty mask.');
end

d0 = zeros(Vmask.dim);
d  = zeros(Vmask.dim);

% if just one contrast is defined, use this contrast and skip interactive selection
if isfield(SPM,'xCon') && numel(SPM.xCon) == 1 && ~isfinite(Ic0)
  Ic0 = 1;
end

% interactively select contrast(s) if necessary
if numel(Ic0)==1 && ~isfinite(Ic0) && numel(SPM.xCon) > 1
  [Ic0,xCon] = spm_conman(SPM,'T&F',Inf,...
    '  Select contrast(s)...',' ',1);
  SPM.xCon = xCon;
end

% go through all contrasts
for con = 1:length(Ic0)
    
  Ic = Ic0(con);
  xCon = SPM.xCon(Ic);
  
  n_perm = job.conspec.n_perm(1);
  if numel(job.conspec.n_perm) > 1
    n_perm_break = job.conspec.n_perm(2);
  end
      
  if ~strcmp(xCon.STAT,'T')
    fprintf('ERROR: Only T-contrasts allowed.\n');
    return
  end
  
  if length(Ic) > 1
    fprintf('ERROR: No conjunction allowed.\n');
    return
  end
  
  fprintf('Use contrast #%d of %s\n',Ic,job.data{1})

  % get contrast and name
  c0 = xCon.c;  
  
  [indi, indj] = find(c0~=0);
  ind_X = unique(indi)';
  xCon.ind_X = ind_X;
  
  % check for contrasts that are defined for columns with subject effects
  if ~isempty(xX.iB)
    if max(ind_X) > min(xX.iB)
      fprintf('ERROR: No contrasts on subjects/block effects allowed.\n');
      return
    end
  end

  % find exchangeability blocks using contrasts without zero values
  exch_blocks   = c0(ind_X,:);
    
  n_exch_blocks = length(ind_X);
    
  % check for exchangeability blocks and design matrix
  if n_exch_blocks == 1
    n_cond = length(find(xX.iH==ind_X)); % check whether the contrast is defined at columns for condition effects
  else
    n_cond = 0;
    n_data_cond = [];
    for k=1:length(xX.iH)
      n_data_cond = [n_data_cond sum(xX.X(:,xX.iH(k)))];
    end
    
    for j=1:n_exch_blocks
      col_exch_blocks = find(c0==exch_blocks(j));
      for k=1:length(col_exch_blocks)
        n_cond = n_cond + length(find(xX.iH==col_exch_blocks(k)));
      end
    end  
    
  end

  % Only regressions designs are currently supported
  if n_cond > 0
    fprintf('ERROR: Only regression designs allowed.\n');
    return
  end

  % No interaction is supported
  if n_exch_blocks > 1
    fprintf('ERROR: No interaction is allowed.\n');
    return
  end

  use_half_permutations = 0;
  % check if sample size is equal for both conditions
  if n_cond == 2
    try
      % repated Anova or F-test don't allow to use only half of the permutions
      if sum(n_data_cond(c0==exch_blocks(1))) == sum(n_data_cond(c0==exch_blocks(2)))
        use_half_permutations = 1;
      else
        use_half_permutations = 0;
      end
    end
  end
  
  ind_exch_blocks = cell(n_exch_blocks,1);
  for j=1:n_exch_blocks
    ind_exch_blocks{j} = find(c0==exch_blocks(j));
  end

  fprintf('\n');
  
  % check design
  switch n_cond
  case 0 % correlation
    label = 1:n_data;
    
    fprintf('Multiple regression design found\n');
  case 1 % one-sample t-test
    fprintf('One sample t-test found\n');
    
    % use exchangeability blocks for labels
    label = zeros(1,n_data);
    for j=1:n_exch_blocks
      for k=1:length(ind_exch_blocks{j})
        label(xX.X(:,ind_exch_blocks{j}(k))~=0) = j;
      end
    end
  otherwise  % Anova with at least 2 groups
    fprintf('Anova found\n');

    % use exchangeability blocks for labels
    label = zeros(1,n_data);
    for j=1:n_exch_blocks
      for k=1:length(ind_exch_blocks{j})
        label(xX.X(:,ind_exch_blocks{j}(k))~=0) = j;
      end
    end
  end

  fprintf('\n')

  % get index for label values > 0
  ind_label = find(label > 0);
  
  n_data_with_contrast = length(ind_label);
  
  % estimate # of permutations
  % Anova/correlation: n_perm = (n1+n2+...+nk)!/(n1!*n2!*...*nk!)
  if n_cond ~=1  % Anova/correlation
    n_perm_full = factorial(n_data_with_contrast);
    single_subject = 0;
    
    for i=1:n_cond
      % check whether only a single subject is in one group
      if length(find(label == i)) == 1
        single_subject = 1;
      end
      n_perm_full = n_perm_full/factorial(length(find(label == i)));
    end
    
    if isnan(n_perm_full)
      % correct number of permutations for large samples when factorial is not working
      if (n_cond == 2) && (single_subject == 1)
        n_perm_full = n_data_with_contrast;
      else
        n_perm_full = realmax;
      end
    end

    % find where data are defined for that contrast
    if ~isempty(find(xX.iH == ind_X(1), 1))
      % first checking whether contrasts are defined for iH
      ind_data_defined = find(any(xX.X(:,xX.iH(ind_X)),2));
    else
      ind_data_defined = find(any(xX.X(:,ind_X),2));
    end
    
    % correct ind_label and n_data_with_contrast using ind_data_defined
    ind_label  = ind_data_defined';
    n_data_with_contrast = length(ind_label);
    
    % and restrict exchangeability block labels to those rows
    exch_block_labels_data_defined = exch_block_labels(ind_data_defined);

    n_perm_full = round(n_perm_full);
    
  else  % one-sample t-test: n_perm = 2^n
    n_perm_full = 2^n_data_with_contrast;
    exch_block_labels_data_defined = exch_block_labels;
    ind_data_defined = ind_label;
  end

  if n_perm_full < n_perm
    fprintf('Warning: Maximum number of possible permutations is lower than defined number of permutations: %d\n',n_perm_full);
  end
  
  n_perm = min([n_perm n_perm_full]);

  fprintf('Number of permutations: %d\n',n_perm);
  
  if use_half_permutations
    fprintf('Equal sample sizes: Use half the number of permutations.\n');
  end
  
  fprintf('Exchangeability block/variable: ');
  fprintf('%d ',unique(cell2mat(ind_exch_blocks)));
  fprintf('\n');
  fprintf('# of conditions: %d\n',n_cond);
           
  X = xX.X(:,ind_X);

  % name of contrast
  c_name0 = deblank(xCon.name);

  c_name = sprintf('%s ',c_name0);

  % obtain new mask that is only defined in regions with significant effects
  P = zeros(VY(1).dim);
  P(ind_mask) = calc_GLM(Y,xX,xCon,ind_mask);
  ind_maskP = P(ind_mask) < 0.001;
  
  Y(~ind_maskP,:) = [];
  
  ind_mask = find(mask > 0 & P < 0.001);
  
  % compute unpermuted d-map
  d0 = calc_NCA(Y,xX,xCon,ind_mask,VY(1).dim);
      
  mask_0   = (d0 == 0);
  mask_1   = (d0 ~= 0);
  mask_P   = (d0 > 0);
  mask_N   = (d0 < 0);
  mask_NaN = (mask == 0);
  found_P  = sum(mask_P(:)) > 0;
  found_N  = sum(mask_N(:)) > 0;
    
  % remove all NaN and Inf's
  d0(isinf(d0) | isnan(d0)) = 0;

  if mesh_detected
    if ~isa(SPM.xVol.G,'gifti')
      % check whether path is correct and file exist
      if ~exist(SPM.xVol.G,'file')
        [pathG,nameG,extG] = spm_fileparts(SPM.xVol.G);
        % use new path
        if ~isempty(strfind(pathG,'_32k'))
          SPM.xVol.G = fullfile(fileparts(which('cat12')),'templates_surfaces_32k',[nameG extG]);
        else
          SPM.xVol.G = fullfile(fileparts(which('cat12')),'templates_surfaces',[nameG extG]);
        end
      end
      SPM.xVol.G = gifti(SPM.xVol.G);
    end
  end

  % prepare output files
  Vt = VY(1);
  Vt.dt(1)    = 16;
  Vt.pinfo(1) = 1;

  %---------------------------------------------------------------
  % save unpermuted t map
  %---------------------------------------------------------------
  name = sprintf('NCAd_%04d',Ic);
  Vt.fname = fullfile(cwd,[name file_ext]);
  Vt.descrip = sprintf('NCAd %04d',Ic);
  Vt = spm_data_hdr_write(Vt);
  spm_data_write(Vt,d0);

  % get largest value
  d0_max    = max(d0(:));
  d0_min    = min(d0(:));
      
  % prepare countings
  dperm       = zeros(size(d));
  d_min       = [];
  d_max       = [];
  d_max_th    = [];

  % general initialization
  try % use try commands to allow batch mode without graphical output
    Fgraph = spm_figure('GetWin','Graphics');
    spm_figure('Clear',Fgraph);
    figure(Fgraph)
  
    h = axes('position',[0.45 0.95 0.1 0.05],'Units','normalized','Parent',...
      Fgraph,'Visible','off');
      
    text(0.5,0.6,c_name,...
      'FontSize',spm('FontSize',10),...
      'FontWeight','Bold',...
      'HorizontalAlignment','Center',...
      'VerticalAlignment','middle')
  
    text(0.5,0.25,spm_str_manip(spm_fileparts(job.data{1}),'a80'),...
      'FontSize',spm('FontSize',8),...
      'HorizontalAlignment','Center',...
      'VerticalAlignment','middle')
  end
  
  % check that label has correct dimension
  sz = size(label);
  if sz(1)>sz(2)
    label = label';
  end
    
  stopStatus = false;
  NCA_progress('Init',n_perm,'Calculating','Permutations');
  
  % update interval for progress bar
  progress_step = max([1 round(n_perm/100)]);

  % Regression design found where contrast is defined for covariate?
  if ~isempty(xX.iC) && all(ismember(ind_X,SPM.xX.iC))
    ind_label_gt0 = find(label(ind_data_defined) > 0);
  else
    ind_label_gt0 = find(label > 0);
  end
  
  unique_labels = unique(label(ind_label_gt0));
  n_unique_labels = length(unique_labels);
  
  perm = 1;
  while perm<=n_perm

    % randomize subject vector
    if perm==1 % first permutation is always unpermuted model
      if n_cond == 1 % one-sample t-test
        rand_label = ones(1,n_data_with_contrast);
        label_matrix = rand_label;
      else % correlation or Anova
        rand_order = ind_label;
        rand_order_sorted = rand_order;
        label_matrix = rand_order;
      end
    else
      % init permutation and
      % check that each permutation is used only once
      if n_cond == 1 % one-sample t-test
        rand_label = sign(randn(1,n_data_with_contrast));
        while any(ismember(label_matrix,rand_label,'rows'))
          rand_label = sign(randn(1,n_data_with_contrast));
        end
      else % correlation or Anova
        
        % permute inside exchangeability blocks only
        rand_order = zeros(1,n_data_with_contrast);
        rand_order_sorted = zeros(1,n_data_with_contrast);
        for k = 1:max(exch_block_labels_data_defined)
          ind_block   = find(exch_block_labels_data_defined == k);
          n_per_block = length(ind_block);
          rand_order(ind_block) = ind_label(ind_block(randperm(n_per_block)));
        end
        
        % go through defined labels and sort inside
        for k=1:n_unique_labels
          ind_block = find(label(ind_label_gt0) == unique_labels(k));
          rand_order_sorted(ind_block) = sort(rand_order(ind_block));
        end

        % check whether this permutation was already used
        count_trials = 0;
        while any(ismember(label_matrix,rand_order_sorted,'rows'))
          count_trials = count_trials + 1;
          
          % stop the permutation loop for too many successless trials for finding 
          % new permutations
          if count_trials > 100000
            fprintf('Stopped after %d permutations because there were too many successless trials for finding new permutations.\n',perm);
            fprintf('Probably there are some missing values for some subjects and the number of maximal permutations was too high.\n');
            n_perm = perm; % stop the permutation loop
            break
          end
          
          for k = 1:max(exch_block_labels_data_defined)
            ind_block   = find(exch_block_labels_data_defined == k);
            n_per_block = length(ind_block);
            rand_order(ind_block) = ind_label(ind_block(randperm(n_per_block)));
          end
          
          % go through defined labels and sort inside
          for k=1:n_unique_labels
            ind_block = find(label(ind_label_gt0) == unique_labels(k));
            rand_order_sorted(ind_block) = sort(rand_order(ind_block));
          end

        end
      end    
    end   
    
    % create permutation set
    Pset = sparse(n_data,n_data);
    if n_cond == 1 % one-sample t-test
      for k=1:n_data_with_contrast
        Pset(ind_label(k),ind_label(k)) = rand_label(k);
      end
    else % correlation or Anova
      for k=1:n_data_with_contrast
        Pset(rand_order_sorted(k),ind_label(k)) = 1;
      end
    end

    % add Stop button after 20 iterations
    try % use try commands to allow batch mode without graphical output
      if perm==21
        hStopButton = uicontrol(Fgraph,...
          'position',[10 10 70 20],...
          'style','toggle',...
          'string','Stop',...
          'backgroundcolor',[1 .5 .5]); % light-red
      end
    
      if perm>=21
        stopStatus = get(hStopButton,'value');
      end
    
      % check Stop status
      if (stopStatus == true)
        fprintf('Stopped after %d iterations.\n',perm);
        break; % stop the permutation loop
      end
    end
      
    % change design matrix according to permutation order
    % only permute columns, where contrast is defined
    Xperm = xX.X;
    Xperm_debug = xX.X;
    Wperm = xX.W;

    Xperm(:,ind_X) = Pset*Xperm(:,ind_X);        
    Xperm_debug(:,ind_X) = Pset*Xperm_debug(:,ind_X);
    
    % correct interaction designs
    % # exch_blocks >1 & # cond == 0 & differential contrast
    if n_exch_blocks >= 2 && n_cond==0 && ~all(exch_blocks(:))
      Xperm2 = Xperm;
      Xperm2(:,ind_X) = 0;
      for j=1:n_exch_blocks
        ind_Xj = find(xX.X(:,ind_X(j)));
        Xperm2(ind_Xj,ind_X(j)) = sum(Xperm(ind_Xj,ind_X),2);
      end
      Xperm = Xperm2;

      Xperm_debug2 = Xperm_debug;
      Xperm_debug2(:,ind_X) = 0;
      for j=1:n_exch_blocks
        ind_Xj = find(xX.X(:,ind_X(j)));
        Xperm_debug2(ind_Xj,ind_X(j)) = sum(Xperm_debug(ind_Xj,ind_X),2);
      end
      Xperm_debug = Xperm_debug2;
    end
    
    if show_permuted_designmatrix
      % scale covariates and nuisance variables to a range 0.8..1
      % to properly display these variables with indicated colors
      if ~isempty(xX.iC)
        val = Xperm_debug(:,xX.iC);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iC) = val;
      end
      if ~isempty(xX.iG)
        val = Xperm_debug(:,xX.iG);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iG) = val;
      end
      if ~isempty(xX.iH) && n_cond==1 % one-sample t-test
        val = Xperm_debug(:,xX.iH);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iH) = val;
      end
      
      % use different colors for indicated columns
      Xperm_debug(:,xX.iH) = 16*Xperm_debug(:,xX.iH);
      Xperm_debug(:,xX.iC) = 24*Xperm_debug(:,xX.iC);
      Xperm_debug(:,xX.iB) = 32*Xperm_debug(:,xX.iB);
      Xperm_debug(:,xX.iG) = 48*Xperm_debug(:,xX.iG);

      if n_cond==1 % one-sample t-test
        for j=1:n_data_with_contrast
          if rand_label(j) > 0
            Xperm_debug(ind_label(j),ind_X) = 60*rand_label(j)*Xperm_debug(ind_label(j),ind_X);
          else
            Xperm_debug(ind_label(j),ind_X) = 56*rand_label(j)*Xperm_debug(ind_label(j),ind_X);
          end
        end
      else % correlation or Anova
        % scale exchangeability blocks also to values 0.8..1
        val = Xperm_debug(:,ind_X);
        ind0 = (val==0);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
      
        % rescue zero entries
        val(ind0) = 0;
      
        Xperm_debug(:,ind_X) = 60*val;
      end

    end
          
    show_plot = 0;
    if use_half_permutations
      if ~rem(perm,progress_step) || ~rem(perm+1,progress_step)
        show_plot = 1;
      end
    else
      if ~rem(perm,progress_step)
        show_plot = 1;
      end
    end

    % display permuted design matrix
    try
      if show_permuted_designmatrix && show_plot
        figure(Fgraph);
        subplot(2,2,3);
        image(Xperm_debug); axis off
        title('Permuted design matrix','FontWeight','bold');
      
        % use different colormap for permuted design matrix
        cmap = jet(64);
      
        % zero values should be always black
        cmap(1,:) = [0 0 0];
        colormap(cmap)
      
        % show legend only once
        if perm <= progress_step
          subplot(2,2,4); axis off
        
          % color-coded legend
          y = 1.0;
          text(-0.2,y, 'Columns of design matrix: ', 'Color',cmap(1, :),'FontWeight','Bold','FontSize',10); y = y - 0.10;
          text(-0.2,y,['Exch. block: ' num2str_short(unique(cell2mat(ind_exch_blocks))')], 'Color',cmap(60,:),'FontWeight','Bold','FontSize',10); y = y - 0.05;
          if ~isempty(xX.iH)
            text(-0.2,y, ['iH - Indicator variable: ' num2str_short(xX.iH)], 'Color',cmap(16,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05; 
          end
          if ~isempty(xX.iC)
            text(-0.2,y, ['iC - Covariate: ' num2str_short(xX.iC)], 'Color',cmap(24,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iB)
            text(-0.2,y, ['iB - Block variable: ' num2str_short(xX.iB)], 'Color',cmap(32,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iG)
            text(-0.2,y, ['iG - Nuisance variable: ' num2str_short(xX.iG)], 'Color',cmap(48,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
        end
      end
    end
    
    % calculate permuted d-map
    if perm == 1
      d    = d0;
    else
      xXperm   = xX;
      xXperm.X = Xperm;        
      xXperm.W = Wperm;

      d = calc_NCA(Y,xXperm,xCon,ind_mask,VY(1).dim);
      
      % remove all NaN and Inf's
      d(isinf(d) | isnan(d)) = 0;
      
    end 
      
    mask_stat_P = mask_1;
    mask_stat_N = mask_1;
    
    % update label_matrix to check for unique permutations
    if use_half_permutations
      if perm>1
        label_matrix = [label_matrix; rand_order_sorted; [rand_order_sorted(label(ind_label) == 2) rand_order_sorted(label(ind_label) == 1)]];
      end
      
      % maximum statistic
      d_max    = [d_max    max(d(mask_stat_P))    -min(d(mask_stat_N))];
      d_min    = [d_min    min(d(mask_stat_N))    -max(d(mask_stat_P))];
      dperm(mask_P)    = dperm(mask_P) + 2*(d(mask_P) >= d0(mask_P));
      dperm(mask_N)    = dperm(mask_N) - 2*(d(mask_N) <= d0(mask_N));
        
    else
      if perm>1
        if n_cond == 1 % one-sample t-test
          label_matrix = [label_matrix; rand_label];
        else
          label_matrix = [label_matrix; rand_order_sorted];
        end
      end

      % maximum statistic
      d_max    = [d_max    max(d(mask_stat_P))];
      d_min    = [d_min    min(d(mask_stat_N))];
      dperm(mask_P)    = dperm(mask_P) + (d(mask_P) >= d0(mask_P));
      dperm(mask_N)    = dperm(mask_N) - (d(mask_N) <= d0(mask_N));
  
    end
      
    % use cummulated sum to find threshold
    st_max    = sort(d_max);

    % find corrected thresholds
    ind_max  = ceil((1-alpha).*length(st_max));
    d_max_th = [d_max_th; st_max(ind_max)];
    if use_half_permutations
      d_max_th = [d_max_th; st_max(ind_max)];
    end
        
    % plot thresholds and histograms      
    try
      if show_plot
        figure(Fgraph);
        axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
        plot_distribution(st_max, d_max_th, 'NCA', alpha, col, 1, d0_max, d0_min);
      end
    end
  
    if numel(job.conspec.n_perm) > 1
      if perm > n_perm_break
        if isempty(find(d0_max > d_max_th(50:end,1), 1))
          fprintf('No FWE-corrected suprathreshold value after %d permutations found\n', n_perm_break);
          perm = n_perm;
        end
      end  
    end
                
    if show_plot
      NCA_progress('Set',perm,Fgraph);
      drawnow
    end
      
    if use_half_permutations  
      perm = perm + 2;
    else
      perm = perm + 1;
    end
  
  end
  
  NCA_progress('Clear',Fgraph);
  
  try
    delete(hStopButton);
    spm_print;
  end
  
  % get correct number of permutations in case that process was stopped
  n_perm = length(d_max);

  %---------------------------------------------------------------
  % corrected threshold based on permutation distribution
  %---------------------------------------------------------------

  % prepare output files
  Vt = VY(1);
  Vt.dt(1)    = 16;
  Vt.pinfo(1) = 1;

  % save ascii file with number of permutations
  name = sprintf('NCAd_%04d',Ic);
  fid = fopen(fullfile(cwd,[name '.txt']),'w');
  fprintf(fid,'%d\n',n_perm);
  fclose(fid);

  %---------------------------------------------------------------
  % save uncorrected p-values for d
  %---------------------------------------------------------------
  name = sprintf('NCAd_log_p_%04d',Ic);
  Vt.fname = fullfile(cwd,[name file_ext]);
  Vt.descrip = sprintf('NCAd %04d',Ic);

  nPtlog10 = zeros(size(d0));

  % estimate p-values
  nPt = dperm/n_perm;
 
  if found_P
    nPtlog10(mask_P) = -log10(nPt(mask_P));
  end
  
  if found_N
    nPt(mask_N) = -nPt(mask_N);
    nPtlog10(mask_N) =  log10(nPt(mask_N));
  end
  
  nPt(mask_0)   = 1;
  nPt(mask_NaN) = NaN;
  nPtlog10(mask_0)   = 0;
  nPtlog10(mask_NaN) = NaN;

  Vt = spm_data_hdr_write(Vt);
  spm_data_write(Vt,nPtlog10);

  %---------------------------------------------------------------
  % save corrected p-values for T
  %---------------------------------------------------------------
  name = sprintf('NCAd_log_pFWE_%04d',Ic);
  Vt.fname = fullfile(cwd,[name file_ext]);
  Vt.descrip = sprintf('NCAd %04d FWE',Ic);

  corrP = zeros(size(d));

  if use_gamma_tail_approximation

    fprintf('Using tail approximation from the Gamma distribution for corrected P-values.\n');

    if tail_approximation_wo_unpermuted_data
      ind_tail = 2:n_perm;
    else
      ind_tail = 1:n_perm;
    end

    if found_P
      [mu,s2,gamm1] = palm_moments(d_max(ind_tail)');
      corrP(mask_P) = palm_gamma(d0(mask_P),mu,s2,gamm1,false,1/n_perm);
    end
    
    if found_N
      [mu,s2,gamm1] = palm_moments(-d_min(ind_tail)');
      corrP(mask_N) = -palm_gamma(-d0(mask_N),mu,s2,gamm1,false,1/n_perm);
    end
    
  else
    if found_P
      for t2 = d_max
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_P) = corrP(mask_P) + (t2 > d0(mask_P)  - tol);
      end
    end
    
    if found_N
      for t2 = d_min
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_N) = corrP(mask_N) - (t2 < d0(mask_N) + tol);
      end
    end
    
    corrP = corrP/n_perm;  
  end

  corrPlog10 = zeros(size(d0));

  if found_P
    corrP(mask_P) = corrP(mask_P);
    corrPlog10(mask_P) = -log10(corrP(mask_P));
  end
  
  if found_N
    corrP(mask_N) = -corrP(mask_N);
    corrPlog10(mask_N) =  log10(corrP(mask_N));
  end

  corrP(mask_0)   = 1;
  corrP(mask_NaN) = NaN;
  corrPlog10(mask_0)   = 0;
  corrPlog10(mask_NaN) = NaN;

  Vt = spm_data_hdr_write(Vt);
  spm_data_write(Vt,corrPlog10);

  %---------------------------------------------------------------
  % save corrected FDR-values for d
  %---------------------------------------------------------------
  name = sprintf('NCAd_log_pFDR_%04d',Ic);
  Vt.fname = fullfile(cwd,[name file_ext]);
  Vt.descrip = sprintf('NCAd %04d FDR',Ic);

  corrPfdr = NaN(size(d));
  corrPfdrlog10 = zeros(size(d0));

  if found_P
    [snP_pos,I_pos] = sort(nPt(mask_P));
    if ~isempty(snP_pos)
      corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
      corrPfdr_pos(I_pos) = corrPfdr_pos;
      corrPfdr(mask_P) = corrPfdr_pos;
      corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
    end
  end

  if found_N
    [snP_neg,I_neg] = sort(nPt(mask_N));
    if ~isempty(snP_neg)
      corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
      corrPfdr_neg(I_neg) = corrPfdr_neg;
      corrPfdr(mask_N) = corrPfdr_neg;
      corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
    end
  end

  corrPfdrlog10(mask_0)   = 0;
  corrPfdrlog10(mask_NaN) = NaN;

  Vt = spm_data_hdr_write(Vt);
  spm_data_write(Vt,corrPfdrlog10);
    
    
end

colormap(gray)

%---------------------------------------------------------------
function plot_distribution(val_max,val_th,name,alpha,col,order,val0_max,val0_min)

corr = 1;

n = length(val_th);
sz_val_max = length(val_max);

% allow other thresholds depending on # of permutations
n_alpha = 3;
if sz_val_max < 1000, n_alpha = 2; end
if sz_val_max <  100, n_alpha = 1; end

% with 20 values we have the lowest possible alpha of 0.05
if sz_val_max >= 20
  alpha  = alpha(1:n_alpha);
  val_th = val_th(:,1:n_alpha);

  [hmax, xmax] = hist(val_max, 100);
      
  subplot(2,2,(2*order)-1)
  
  h = bar(xmax,hmax);
  set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);

  max_h = max(hmax);
  lim_x = xlim;

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color',[0.3333 1 0],'LineWidth',2);
  text(0.95*val0_max,0.95*max_h,'Max. observed value ',...
    'Color',[0.3333 1 0],'HorizontalAlignment','Right','FontSize',8)
    
  % plot sign-flipped minimum observed value for unpermuted model
  if val0_min < 0
    hl = line([-val0_min -val0_min], [0 max_h]);
    set(hl,'Color',[0 0.6667 1],'LineWidth',2);
    text(0.95*val0_max,0.85*max_h,'Max. observed value (inverse contrast) ',...
      'Color',[0 0.6667 1],'HorizontalAlignment','Right','FontSize',8)
  end
  
  % plot thresholds
  for j=1:n_alpha
    hl = line([val_th(n,j) val_th(n,j)], [0 max_h]);
    set(hl,'Color',col(j,:),'LineStyle','--');
    text(0.95*lim_x(2),(0.4+0.1*j)*max_h,['p<' num2str(alpha(j))],...
      'Color',col(j,:),'HorizontalAlignment','Right','FontSize',8)
  end
  
  ylabel('Frequency');
  xlabel(['Max ' name]);
  if corr
    title(['Distribution of maximum ' name],'FontWeight','bold');
  else
    title(['Distribution of ' name],'FontWeight','bold');
  end
  
  subplot(2,2,2*order)
  
  val_min = min(min(val_th(1:n,:)));
  val_max = max(max(val_th(1:n,:)));
  if val_max/val_min > 10
    hp = semilogy(1:n,val_th(1:n,:));
    yl = log10(ylim);
    ylim(10.^[floor(yl(1)) ceil(yl(2))])
  else
    hp = plot(1:n,val_th(1:n,:));
  end
  
  for j=1:n_alpha
    set(hp(j),'Color',col(j,:));
  end
  
  if corr
    title(['Corrected threshold of ' name],'FontWeight','bold')
  else
    title(['Uncorr. threshold of ' name],'FontWeight','bold')
  end
  ylabel('Threshold')
  xlabel('Permutations')   
end

%---------------------------------------------------------------
function d = calc_NCA(Y,xX,xCon,ind_mask,dim)
% compute effect size of the necessary condition
%
% Y        - masked data as vector
% xX       - design structure
% xCon     - contrast structure
% ind_mask - mask image
% dim      - image dimension
%
% Output:
% d        - effect size d

c = xCon.c;
ind_c = c ~= 0;
if sum(ind_c) > 1
  error('Only regression designs supported.');
end

X = xX.W*xX.X(:,ind_c);

d = zeros(dim);
d(ind_mask) = NCA(X,Y');

if c(ind_c) < 0
  d = -d;
end

%---------------------------------------------------------------
function P = calc_GLM(Y,xX,xCon,ind_mask)
% compute T- or F-statistic using GLM
%
% Y        - masked data as vector
% xX       - design structure
% xCon     - contrast structure
% ind_mask - index of mask image
%
% Output:
% T        - P-values

c = xCon.c;

X = xX.W*xX.X;
pKX = pinv(X);

n_data = size(X,1);

Beta = Y*pKX';

res0 = Beta*(single(X'));
res0 = Y - res0; %-Residuals
res0 = res0.^2;
ResSS = double(sum(res0,2));
clear res0

trRV = n_data - rank(xX.X);
ResMS = ResSS/trRV;
%-Modify ResMS (a form of shrinkage) to avoid problems of very low variance
ResMS  = ResMS + 1e-3 * max(ResMS(isfinite(ResMS)));

if strcmp(xCon.STAT,'T')
  Bcov = pKX*pKX';
  con = Beta*c;

  T = con./(eps+sqrt(ResMS*(c'*Bcov*c)));
  P = 1-spm_Tcdf(T,trRV);
else
  error('F-test not supported');
end


%---------------------------------------------------------------
function xX = correct_xX(xX)
% sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i, 1))
    Xsum(:,i) = sum(X(:,1:i),2);
  end
end

% find columns where all entries are constant except zeros entries
% that indicate columns with covariates and nuisance variables
ind = find(any(diff(Xsum))==0 & sum(Xsum)>0);

% no more than 2 factors expected
if length(ind) > 2
  error('Weird design was found that cannot be analyzed correctly.');
end

% correction is only necessary if 2 factors (iH/iB) were found
if length(ind) > 1
  iF = cell(length(ind),1);

  j = 1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j),  j = j + 1; end

  for i=j:length(ind)
    iF{i} = [j:ind(i)];
  
    j = ind(i)+1;
    % skip columns with covariates and nuisance variables
    while find(iCG==j), j = j + 1; end
  end
  
  % not sure whether this will always work but usually iB (subject effects) should be larger than iH (time effects)
%  if length(iF{1}) > length(iF{2})
if 0 % will be probably not always correct 
    xX.iB = iF{1};
    xX.iH = iF{2};
  else
    xX.iB = iF{2};
    xX.iH = iF{1};
  end
end

%---------------------------------------------------------------
function str = num2str_short(num)
% get shorther strings for continuous numbers with length > 4

if length(num) > 4
  % check whether vector consist of continuous numbers
  if all(diff(num)==1)
    str = [num2str(num(1)) ':' num2str(num(end))];
  else
    str = num2str(num);
  end
else
  str = num2str(num);
end

%---------------------------------------------------------------
function varargout = palm_moments(varargin)
% For a statistic that can be expressed as trace(A*W), for
% a sample size of n observations, this function returns the
% expected first three moments of the permutation distribution,
% without actually computing any permutation.
%
% [mu,sigsq,gamm1,gamm2] = palm_moments(G)
% [mu,sigsq,gamm1]       = palm_moments(A,W,n)
%
% Inputs:
% - G : A PxV array of observations of the G random variable.
%       The moments are unbiased and run along the 1st dimension.
%       Typical case is P = number of permutations and V = 
%       number of tests, e.g., voxels.
% - A : A square matrix for a multivariate proper statistic,
%       or a vector of A values for various univariate tests.
% - W : A square matrix for a multivariate proper statistic,
%       or a vector of W values for various univariate tests.
% - n : Sample size on which the statistic is based.
%
% Outputs:
% - mu    : Sample mean.
% - sigsq : Sample variance (unbiased).
% - gamm1 : Sample skewness (unbiased).
% - gamm2 : Sample kurtosis (unbiased).
%
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% For the estimators using trace(AW), the references are:
% * Kazi-Aoual F, Hitier S, Sabatier R, Lebreton J-D. Refined
%   approximations to permutation tests for multivariate
%   inference. Comput Stat Data Anal. 1995;20(94):643-656.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;4:497-511.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin == 1
    
    % For a set of values of the random variable G, return the
    % first 4 moments.
    
    % Mean
    G     = varargin{1};
    n     = size(G,1);
    mu    = sum(G,1)/n;
    
    % Variance
    G0    = bsxfun(@minus,G,mu);
    ssq   = sum(G0.^2,1);
    sigsq = (ssq/(n-1));
    
    % Skewness
    s2    = ssq/n; % biased variance
    m3    = sum(G0.^3,1)/n;
    gamm1 = m3./s2.^1.5;
    gamm1 = gamm1 * sqrt((n-1)/n)*n/(n-2); % unbiased skewness

    % Kurtosis (normal dist = 3)
    if nargout == 4
        m4    = sum(G0.^4,1)/n;
        gamm2 = (m4./s2.^2);
        gamm2 = ((n+1)* gamm2 -3*(n-1))*(n-1)/((n-2)*(n-3))+3; % unbiased kurtosis
    else
        gamm2 = [];
    end
    
elseif nargin == 3
    
    % Compute the first three moments of the permutation distribution of
    % the statistic G = trace(AW), for n subjects, using the method in
    % Kazi-Aoual et al (1995). The variable names follow roughly the same
    % as in the paper.
    
    % Take inputs
    A = varargin{1};
    W = varargin{2};
    n = varargin{3};
    
    % If A and W are truly multivariate (i.e., square matrices), do as in
    % the original paper. Otherwise, make simplifications as these are all
    % scalars.
    if size(A,1) == size(A,2)
        
        % Some auxiliary variables for ET2:
        T    = trace(A);
        T_2  = trace(A^2);
        S_2  = sum(diag(A).^2);
        Ts   = trace(W);
        Ts_2 = trace(W^2);
        Ss_2 = sum(diag(W).^2);
        
        % Some auxiliary variables for ET3:
        T_3  = trace(A^3);
        S_3  = sum(diag(A).^3);
        U    = sum(A(:).^3);
        R    = diag(A)'*diag(A^2);
        B    = diag(A)'*A*diag(A);
        Ts_3 = trace(W^3);
        Ss_3 = sum(diag(W).^3);
        Us   = sum(W(:).^3);
        Rs   = diag(W)'*diag(W^2);
        Bs   = diag(W)'*W*diag(W);
        
    else
        
        % Some auxiliary variables for ET2:
        T    = A;
        T_2  = A.^2;
        S_2  = T_2;
        Ts   = W;
        Ts_2 = W.^2;
        Ss_2 = Ts_2;
        
        % Some auxiliary variables for ET3:
        T_3  = A.^3;
        S_3  = T_3;
        U    = T_3;
        R    = T_3;
        B    = T_3;
        Ts_3 = W.^3;
        Ss_3 = Ts_3;
        Us   = Ts_3;
        Rs   = Ts_3;
        Bs   = Ts_3;
    end
    
    % E(T):
    mu = T.*Ts/(n-1);
    
    % V(T):
    sigsq = 2*((n-1)*T_2-T.^2).*((n-1)*Ts_2-Ts.^2) / (n-1)^2/(n+1)/(n-2) ...
        + (n*(n+1)*S_2-(n-1)*(T.^2+2*T_2)) .* (n*(n+1)*Ss_2-(n-1)*(Ts.^2+2*Ts_2)) ...
        / (n+1)/n/(n-1)/(n-2)/(n-3);
    
    % E(T^3):
    ET3 = ...
        n^2*(n+1)*(n^2+15*n-4)*S_3.*Ss_3 ...
        + 4*(n^4-8*n^3+19*n^2-4*n-16)*U.*Us ...
        + 24*(n^2-n-4)*(U.*Bs+B.*Us) ...
        + 6*(n^4-8*n^3+21*n^2-6*n-24)*B.*Bs ...
        + 12*(n^4-n^3-8*n^2+36*n-48)*R.*Rs ...
        + 12*(n^3-2*n^2+9*n-12)*(T.*S_2.*Rs + R.*Ts.*Ss_2) ...
        + 3*(n^4-4*n^3-2*n^2+9*n-12)*T.*Ts.*S_2.*Ss_2 ...
        + 24*( (n^3-3*n^2-2*n+8)*(R.*Us+U.*Rs) ...
        + (n^3-2*n^2-3*n+12)*(R.*Bs+B.*Rs) ) ...
        + 12*(n^2-n+4)*(T.*S_2.*Us+U.*Ts.*Ss_2) ...
        + 6*(2*n^3-7*n^2-3*n+12)*(T.*S_2.*Bs+B.*Ts.*Ss_2) ...
        - 2*n*(n-1)*(n^2-n+4)*( (2*U+3*B).*Ss_3+(2*Us+3*Bs).*S_3 ) ...
        - 3*n*(n-1)^2*(n+4)*( (T.*S_2+4*R).*Ss_3+(Ts.*Ss_2+4*Rs).*S_3 ) ...
        + 2*n*(n-1)*(n-2)*( (T.^3+6*T.*T_2+8*T_3).*Ss_3 ...
        + (Ts.^3+6*Ts.*Ts_2+8*Ts_3).*S_3 ) ...
        + T.^3.*((n^3-9*n^2+23*n-14)*Ts.^3+6*(n-4).*Ts.*Ts_2+8*Ts_3) ...
        + 6*T.*T_2.*((n-4)*Ts.^3+(n^3-9*n^2+24*n-14)*Ts.*Ts_2+4*(n-3)*Ts_3) ...
        + 8*T_3.*(Ts.^3+3*(n-3).*Ts.*Ts_2+(n^3-9*n^2+26*n-22)*Ts_3) ...
        - 16*(T.^3.*Us+U.*Ts.^3)-6*(T.*T_2.*Us+U.*Ts.*Ts_2)*(2*n^2-10*n+16) ...
        - 8*(T_3.*Us+U.*Ts_3)*(3*n^2-15*n+16)-(T.^3.*Bs+B.*Ts.^3) ...
        * (6*n^2-30*n+24)-6*(T.*T_2.*Bs+B.*Ts.*Ts_2)*(4*n^2-20*n+24) ...
        - 8*(T_3.*Bs + B.*Ts_3)*(3*n^2-15*n+24) ...
        - (n-2)*( 24*(T.^3.*Rs+R.*Ts.^3)+6*(T.*T_2.*Rs+R.*Ts.*Ts_2)*(2*n^2-10*n+24) ...
        + 8*(T_3.*Rs+R.*Ts_3)*(3*n^2-15*n+24)+(3*n^2-15*n+6) ...
        .* (T.^3.*Ts.*Ss_2+T.*S_2.*Ts.^3) ...
        + 6*(T.*T_2.*Ts.*Ss_2+T.*S_2.*Ts.*Ts_2)*(n^2-5*n+6) ...
        + 48*(T_3.*Ts.*Ss_2+T.*S_2.*Ts_3) );
    ET3 = ET3/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5);
    
    % The coefficient "3" below is missing from Kazi-Aoual et al (1995), but it
    % is shown in the Supplementary Information of Minas and Montana (2014).
    gamm1 = (ET3 - 3*mu.*sigsq - mu.^3)./sigsq.^1.5;
    gamm2 = [];
else
    error('Incorrect number of arguments');
end

% Return results
varargout{1} = mu;
varargout{2} = sigsq;
varargout{3} = gamm1;
varargout{4} = gamm2;

%---------------------------------------------------------------
function pvals = palm_gamma(G,mu,sigsq,gamm1,rev,prepl)
% Return the p-values for a Gamma distribution, parameterised by
% its first three moments.
%
% pvals = palm_gamma(G,mu,s2,gamm1,rev)
% 
% Inputs:
% - G     : Statistics for which p-values are to be computed.
% - mu    : Distribution mean.
% - sigsq : Distribution standard deviation.
% - gamm1 : Distribution skewness.
% - rev   : Use if lower values of the statistic are evidence in
%           favour of the alternative.
% - prepl : Replacement for what otherwise would be zero p-values
%           in case of poor fits (e.g., statistic falls into the
%           part of the distribution that has pdf=0. In these cases
%           the p-value can be 1 or 1/(#perm) depending on which
%           tail and the sign of the skewness.
%
% Outputs:
% - pvals : p-values.
% 
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% Other references:
% * Mielke PW, Berry KJ, Brier GW. Application of Multi-Response
%   Permutation Procedures for Examining Seasonal Changes in
%   Monthly Mean Sea-Level Pressure Patterns. Mon Weather Rev.
%   1981;109(1):120-126.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;7(6):450-470.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Note that there are no argument checking for speed, but
% sizes of all inputs need to be the same, or the moments need to
% be all scalars.

if gamm1 == 0
    
    % If not skewed, use a normal approximation.
    G     = (G - mu)./sigsq.^.5;
    pvals = erfc(G/sqrt(2))/2;
    
else
    
    % Standardise G, so that all becomes a function of the skewness.
    G    = (G - mu)./sigsq.^.5;
    
    % Gamma distribution parameters (Minas & Montana, 2014).
    kpar = 4/gamm1.^2;
    tpar = gamm1/2;
    cpar = -2/gamm1;
     
    % Actual p-value. If there are negatives here, the probability can
    % have an imaginary part, which is dealt with later.
    if rev
        if gamm1 > 0
            tail = 'lower';
        else
            tail = 'upper';
        end
    else
        if gamm1 > 0
            tail = 'upper';
        else
            tail = 'lower';
        end
    end
    pvals = gammainc((G-cpar)./tpar,kpar,tail);
    
    % Deal with imaginary parts.
    if ~ isreal(pvals)
        iidx = imag(pvals) ~= 0;
        if rev
            if gamm1 > 0
                pvals(iidx) = prepl;
            else
                pvals(iidx) = 1;
            end
        else
            if gamm1 > 0
                pvals(iidx) = 1;
            else
                pvals(iidx) = prepl;
            end
        end
    end
end


