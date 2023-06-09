function NCA_estimate = tbx_cfg_NCA
% SPM Configuration file for NCA estimate
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________


addpath(fileparts(which(mfilename)));

% try to estimate number of processor cores
try
  if strcmpi(spm_check_version,'octave')
    numcores = nproc;
  else
    numcores = feature('numcores');
  end

  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
% or for Octave
if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end

%_______________________________________________________________________
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.hidden  = numcores <= 1 || isdeployed;
nproc.help    = {
    'In order to use multi-threading the NCA job with multiple SPM.mat files can be split into separate processes that run in the background. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process might need a large amount of RAM, which should be considered to choose the appropriate number of processes.'
    ''
    'Please further note that additional modules in the batch can now be used because the processes are checked every minute.'
  };

% ---------------------------------------------------------------------
% data Select SPM.mat
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Select SPM.mat';
data.help    = {'Select the SPM.mat files that contain the design specification from a previous (parametric) estimation, where all required contrasts are already specified.'};
data.filter  = 'mat';
data.ufilter = '^SPM\.mat$';
data.num     = [1 Inf];

% ---------------------------------------------------------------------
% mask Select mask to restrict analysis
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Select additional mask';
mask.help    = {'Select an additional mask image or surface to restrict your analysis. As default the mask in the analysis folder is used. Here you can select a mask to additionally restrict the analysis to regions of interest (i.e. small volume/surface correction).'};
if strcmp(spm('ver'),'SPM12')
  mask.filter  = {'image','mesh'};
else
  mask.filter  = {'image'};
end
mask.val     = {''};
mask.ufilter = '.*';
mask.num     = [0 1];

% ---------------------------------------------------------------------
% titlestr Results Title
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results title';
titlestr.help    = {'Heading on results page - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];

% ---------------------------------------------------------------------
% contrasts Contrast
% ---------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast index';
contrasts.help    = {'Index(es) of contrast according to the contrast manager.'
                     ''
                     'Each contrast in SPM is indicated by a sequential number that is displayed in the first column of the contrast manager.'
                     ''
                     'You can enter one or more contrasts. If only one number is entered, and this number is "Inf", you can select one or more contrasts interactively using the contrast manager.'
                     ''
                     'Do not define here the contrast itself. This should be done in the contrast manager, that is automatically called if "Inf" is kept as entry.'
}';
contrasts.strtype = 'e';
contrasts.val     = {Inf};
contrasts.num     = [1 Inf];

% ---------------------------------------------------------------------
% number of permutations
% ---------------------------------------------------------------------
n_perm         = cfg_entry;
n_perm.tag     = 'n_perm';
n_perm.name    = 'Number of permutations';
n_perm.help    = {'In order to obtain reliable estimates you need about 5000-10000 permutations.'
                  ''
                  'There is also an option to interrrupt the permutation process and to save the results at this step to take a first look on your results.'
                  ''
                  'If the number of maximal possible permutations is smaller, then this number is used resulting in an exact permutation test.'
                  ''
                  'Please note, that a tail approximation is finally used to estimate the corrected p-values. Thus, there is no dependency anymore between the lowest achievable p-value and the number of permutations as in previous versions.'
}';
n_perm.strtype = 'e';
n_perm.val     = {5000};
n_perm.num     = [1 Inf];

% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts n_perm};
conspec.help    = {''};

% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
NCA_estimate          = cfg_exbranch;
NCA_estimate.tag      = 'NCA_estimate';
NCA_estimate.name     = 'Estimate NCA';
NCA_estimate.val      = {data nproc mask conspec};
NCA_estimate.help     = {''};
NCA_estimate.prog     = @NCA_estimate_stat;
