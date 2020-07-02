function percentCorrect = svm_compare_scenes(theOI, theMosaic, scene1, scene2, add_poisson_noise)

    tic;
    %% Compute the optical images to the two scenes
    theOIscene1 = oiCompute(theOI, scene1);
    theOIscene2 = oiCompute(theOI, scene2);
    disp('Time: Computed 2 optical images')
    toc;
                    
    tic;
    %% Compute the mean cone excitation responses to the two stimuli
    coneExcitationsCond1 = theMosaic.compute(theOIscene1);
    coneExcitationsCond2 = theMosaic.compute(theOIscene2);
    %saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2); 
    disp('Time: Computed 2 excitations')
    toc;
    
    randomSeed = rng;
    %saveparfor_rseed(savename_coneresp, randomSeed); 
    
    tic;
    if add_poisson_noise
        [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(coneExcitationsCond1, coneExcitationsCond2);
        %saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances);
    end
    disp('Time: Added poisson noise')
    toc;
    
    tic;
    percentCorrect = svm_pca(theMosaic, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances);
    disp('Time: Fit excitations with PCA-SVM model (1024 instances per alternative, 60 dimension PCA)')
    toc;
    %saveparfor_svm(savename_coneresp, SVMpercentCorrect);
end

function load_previous_excitations(savename_coneresp,nSVMrep)
    [runsvm, nrunadd] = checkforsvm(savename_coneresp, nSVMrep);
    if runsvm
        if isfile(savename_coneresp_instances)
            stemp_instances = load(savename_coneresp_instances);
            coneExcitationCond1_noisyInstances = stemp_instances.coneExcitationCond1_noisyInstances;
            coneExcitationCond2_noisyInstances = stemp_instances.coneExcitationCond2_noisyInstances;
        else
            stemp = load(savename_coneresp);
            [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(stemp.coneExcitationsCond1, stemp.coneExcitationsCond2);
            saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances);
        end
        SVMpercentCorrect = [];
        for rep = 1:nrunadd
            SVMpercentCorrect = [SVMpercentCorrect, svm_pca(theMosaic, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances)];
        end

        saveparfor_svm(savename_coneresp, SVMpercentCorrect);
    end
end

function [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(coneExcitationsCond1, coneExcitationsCond2)

    % Compute some noisy instances of cone mosaic excitations
    nInstances = 1024;
    [nr1, nc1] = size(coneExcitationsCond1);
    [nr2, nc2] = size(coneExcitationsCond2);
    coneExcitationCond1_noisyInstances = nan(nInstances, nr1, nc1); 
    coneExcitationCond2_noisyInstances = nan(nInstances, nr2, nc2); 
    for i = 1:nInstances
        coneExcitationCond1_noisyInstances(i,:,:) = poissrnd(coneExcitationsCond1);
        coneExcitationCond2_noisyInstances(i,:,:) = poissrnd(coneExcitationsCond2);
    end
end

function saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2)

fprintf('Saving %s. \n', savename_coneresp);
save(savename_coneresp, 'coneExcitationsCond1', 'coneExcitationsCond2');

end

function saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances)
    save(savename_coneresp_instances, 'coneExcitationCond1_noisyInstances', 'coneExcitationCond2_noisyInstances', '-v7.3');
end

function saveparfor_rseed(savename_coneresp, randomSeed)
    save(savename_coneresp, 'randomSeed', '-append');
end

function saveparfor_svm(savename_coneresp, SVMpercentCorrect)
    fprintf('%s: %.2f \n', savename_coneresp, mean(SVMpercentCorrect));
    if ismember(who('-file', savename_coneresp), 'SVMpercentCorrect')
        SVMpercentCorrect_prev = load(savename_coneresp, 'SVMpercentCorrect');
        SVMpercentCorrect = [SVMpercentCorrect_prev, SVMpercentCorrect];
    end
    save(savename_coneresp, 'SVMpercentCorrect', '-append');
end

function [runsvm, nrunadd] = checkforsvm(savename_coneresp, nSVMrep)
    if ismember(who('-file', savename_coneresp), 'SVMpercentCorrect')
        SVMpercentCorrect = load(savename_coneresp, 'SVMpercentCorrect'); 
    end 
    if length(SVMpercentCorrect) < nSVMrep
        runsvm = true;
        nrunadd = nSVM-length(SVMpercentCorrect);
        if ~ismember(who('-file', savename_coneresp), 'randomSeed')
            randomSeed = rng;
            saveparfor_rseed(savename_coneresp, randomSeed);
        end
        fprintf('%d SVM run exists. Adding %d more runs... \n', length(SVMpercentCorrect), nrunadd);
    else
        runsvm = false;
        nrunadd = 0;
        fprintf('%d SVM run exists - this simulation is complete. \n', length(SVMpercentCorrect)); 
    end 
end 
