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


%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
%-------------------------------------------------%
resultdir = 'Results';
foldername = 'ConeExcitationInstances_SPDcorrectedScaled_60PCA_1024Instances';
KLMSdensity = {[0 13/14 1/14 0]', [0 5/6 1/6 0]', [0 2.8/3.8 1/3.8 0]', [0 1.8/2.8 1/2.8 0]', [0.0 0.5 0.5 0.0]', [0 1/2.8 1.8/2.8 0]', [0 1/3.8 2.8/3.8, 0]', [0 1/6 5/6 0]', [0 1/14 13/14 0]'};
coltype_set = {[0 0], [1 1], [0 1]};     % isochromatic, isoluminant, isochromatic vs. isoluminant
ort_set = {[0 1], [0 1], [1 1]};         % diff orts, diff orts, same ort
sf_set = [4, 8, 16, 24, 32, 48, 64]; 
nSF = length(sf_set); %[10 30 60]; %linspace(5, 50, nSF);             % for each of the three above
nContrast = 15; contrast_set = 10.^linspace(log10(0.001), log10(0.08), nContrast); %10.^linspace(log10(0.03), log10(0.1), nContrast); %linspace(0.001, 0.05, nContrast); % for each of the three above
nSVMrep = 1; 
%-------------------------------------------------%


% Making dir to save all the simulation results
if ~isfolder(resultdir)
    mkdir(resultdir);
end

parpool('local',2); % More parallel workers often leads to memory problems

for mos = 1:length(KLMSdensity)
    
    this_KLMSdensity = KLMSdensity{mos};
    
    % Making dir to save cone excitation instances
    mosaicdir = fullfile(resultdir, ['mosaicCond', num2str(mos)]);
    if ~isfolder(mosaicdir)
        mkdir(mosaicdir);
    end
    
    savename_mosaic = fullfile(mosaicdir, ['mosaic_L', num2str(this_KLMSdensity(2)*10), ...,
                                                  'M', num2str(this_KLMSdensity(3)*10), ...,
                                                  'S', num2str(this_KLMSdensity(4)*10), '.mat']);
    
    if isfile(savename_mosaic)
        load(savename_mosaic);
        fprintf('Loading file: %s \n', savename_mosaic);
        if ~strcmp(theMosaic.noiseFlag, 'none')
            fprintf('The noiseFlag of this mosaic is set to %s. Turning off the noise. /n', theMosaic.noiseFlag);
            theMosaic.noiseFlag = 'none';
        end
    else
        theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
            'fovDegs', sceneFov, ...                    % match mosaic width to stimulus size
            'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
            'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
            'integrationTime', 10/1000, ...             % 30 msec integration time
            'maxGridAdjustmentIterations', 50, ...
            'spatialDensity', this_KLMSdensity, ...     % terminate iterative lattice adjustment after 50 iterations
            'noiseFlag', 'none');                       % For efficiency, we add Poisson noise after-the-fact. Only works
                                                        % for no eye movements, etc.
        save(savename_mosaic, 'theMosaic', 'this_KLMSdensity');
    end
    
    % Making dir to save cone excitation instances
    conerespdir = fullfile(mosaicdir, foldername);
    if ~isfolder(conerespdir)
        mkdir(conerespdir);
        mkdir(fullfile(conerespdir, 'noisyInstances'));
    end
    
    condIdPerColtype = [floor((0:nSF*nContrast-1)/nContrast)' + 1, mod(0:nSF*nContrast-1, nContrast)' + 1];

    for oi = 1:2
        % For oi==1, do LCA off, everything OFF
        % For oi==2, do LCA on, normal human wvf
        if oi == 2
            theOI = theOI_with_lca;
            disp('Now using the normal human WVF OI.');
        end
        
        for exp = 1:length(coltype_set) 
            
            this_coltype = coltype_set{exp};
            this_ort     = ort_set{exp};
            fprintf('Experimental condition %d. \n', exp);
            
            parfor cnd = 1:size(condIdPerColtype, 1)
                close all;
                
                this_sf       = sf_set(condIdPerColtype(cnd, 1));
                this_contrast = contrast_set(condIdPerColtype(cnd, 2));
                
                savename_coneresp = fullfile(conerespdir, ['coneExcitation_noiseOff_oi', num2str(oi), '_exp', num2str(exp), '_SF_', num2str(this_sf), '_contr_', num2str(this_contrast), '.mat']);
                [path, fname, ext] = fileparts(savename_coneresp);
                savename_coneresp_instances = fullfile(path, 'noisyInstances', ['noisyInstances_', fname, ext]);
                
                if isfile(savename_coneresp)
                    fprintf('This cone excitation instance already exists. Loading saved file... \n');
                    load_previous_excitations(savename_coneresp,nSVMrep);
                else
                    %% ------ LOOP FOR EACH CONDITION FROM HERE ------ %%
                    scene1 = generateGaborSceneAO(presentationDisplay, this_coltype(1), this_ort(1), this_sf, this_contrast); % inputs: (display, coltype, ort, sf, contrast)
                    scene2 = generateGaborSceneAO(presentationDisplay, this_coltype(2), this_ort(2), this_sf, this_contrast);
                    %% Compute the mean cone excitation responses to the stimulus

                    svm_compare_scenes(theOI, theMosaic, scene1, scene2, 1); % Last param: add_poisson_noise (Function should add noise)
                end
            end
        end
    end
end
