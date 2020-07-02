clear all;
close all;

rng;

[theOI_with_lca,theOI]=make_optics();

%% Create the scene
presentationDisplay = displayCreate('AOSim-Seattle_SPDcorrected_Scaled');

scene_sample = generateGaborSceneAO(presentationDisplay, 1, 1, 1, 1); % just to get the fov for mosaic generation
sceneFov = 1.1;%sceneGet(scene_sample, 'fov');
% ok, if this value is smaller than that size of the scene, the mosaic
% generation gets error. Here as a quick remedy artificially giving a
% slightly higer fov.


resultdir = 'Results';
foldername = 'ConeExcitationInstances_SPDcorrectedScaled_60PCA_1024Instances';
KLMSdensity1 = [0 2/3 1/3 0]'; % L:M:S = 2:1:0

% Specify 2AFC pair for discrimination task:
colors = [1 1];     		% Both are Red/Green counterphase gratings
orientations = [0 1];  		% Horizontal vs. Vertical

sf = 16;
nContrast = 3; contrast_set = [0.001,0.0125,0.05]; 
results = [];

conerespdir=resultdir;
if ~isfolder(resultdir)
    mkdir(resultdir);
end

tic;
theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
    'fovDegs', sceneFov, ...                    % match mosaic width to stimulus size
    'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
    'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
    'integrationTime', 10/1000, ...             % 30 msec integration time
    'maxGridAdjustmentIterations', 50, ...
    'spatialDensity', KLMSdensity1, ...     % terminate iterative lattice adjustment after 50 iterations
    'noiseFlag', 'none');                       % For efficiency, we add Poisson noise after-the-fact. Only works
						% for no eye movements, etc.
disp('Time: Generate single mosaic, 1deg')
toc;
            
for ncon = 1:nContrast
	close all;
	
	this_contrast = contrast_set(ncon);
	
	savename_coneresp = fullfile(conerespdir, ['coneExcitation_noiseOff_oi', num2str(1), '_exp', num2str(1), '_SF_', num2str(sf), '_contr_', num2str(this_contrast), '.mat']);
	[path, fname, ext] = fileparts(savename_coneresp);
	savename_coneresp_instances = fullfile(path, 'noisyInstances', ['noisyInstances_', fname, ext]);
	
	tic;
	scene1 = generateGaborSceneAO(presentationDisplay, colors(1), orientations(1), sf, this_contrast); 
	scene2 = generateGaborSceneAO(presentationDisplay, colors(2), orientations(2), sf, this_contrast);
    	disp('Time: Generate 2 Gabor Scenes')
	toc;

	pct_correct=svm_compare_scenes(theOI, theMosaic, scene1, scene2, 1); % Last param: add_poisson_noise (Function should add noise)
	results = [results,pct_correct];
end

% For contrasts of [0.001,0.0125,0.05], results are approx. 50%, 80%, 100%
disp('Proportions correct for 3 contrasts:')
disp(results);
