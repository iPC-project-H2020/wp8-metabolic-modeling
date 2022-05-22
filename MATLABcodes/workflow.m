cd ' ' %%% Set the working directory here
[num,txt,raw] = xlsread('job_profiles.xlsx');

njobs = size(raw);
varNames = raw(1,:);
raw = cell2table(raw, 'VariableNames', varNames);

for ii = 2:njobs(1)
    jobInfo = raw(ii,:);
    subdir = constructDirectory(jobInfo);
    subdir = subdir{1};
    %%% pre-set result directories
    FBAdir = strcat(pwd, '\FBA_solution', subdir);
    ValToReact_dir = strcat(pwd, '\ValToReact', subdir);
    csGEMdir = strcat(pwd, '\csGEM', subdir);
    if not(isfolder(FBAdir))
        mkdir(FBAdir);
        mkdir(csGEMdir);
        mkdir(ValToReact_dir);
    else 
        continue
    end
    %%% set and select cluster profile data
    clusterProfile_dir = strcat(pwd, '\cluster_profiles', subdir);
    dsClusterProfiles = datastore(clusterProfile_dir,'Type', 'tabulartext');
    dsProfilesTransformed = transform(dsClusterProfiles, @addFilenameToData, ...
    'IncludeInfo', true);
    %%% select, load and preprocess GEM model
    GEMdir = strcat(pwd, '\GEM\', jobInfo.model_name{1});
    modelFileName = strcat(GEMdir, '\', jobInfo.model_name{1}, '.mat');
    load(modelFileName);
    [Nmets, ~ ]= size(model.mets);
    model.osenseStr = 'max';
    
    while hasdata(dsClusterProfiles)
        T = table(model.rxns);
        ClusterProfile = read(dsClusterProfiles);
        fileName = read(dsProfilesTransformed);
        fileName = fileName.Filename{1}
        fileName = split(fileName, '.txt');
        fileName = fileName{1};
    
        exprData.value = ClusterProfile.Var2(find(ClusterProfile.Var2~=0));
        exprData.gene = ClusterProfile.Var1(find(ClusterProfile.Var2~=0));
        [gene_id, gene_expr] = findUsedGenesLevels(model, exprData);
        [expressionRxns parsedGPR] = mapExpressionToReactions(model, exprData);
    
        expReact = table(expressionRxns);
        expReact.Properties.VariableNames = {fileName};
        T = [T, expReact];
        writetable(T, strcat(ValToReact_dir, fileName, '.txt')) 
    
        thresHold = quantile(expressionRxns, 0.25) ;
        changeCobraSolver('glpk');
    
        specModel = GIMME(model, expressionRxns, thresHold, 0.9);
        ctName = strcat(csGEMdir, 'csGEM_' ,fileName, '.mat');
        save([ctName], 'specModel');    
        SolBiomass = optimizeCbModel(specModel); 
        FBAsolBiomass = specModel.rxns;

        writetable([FBAsolBiomass,array2table(SolBiomass.x)], strcat(FBAdir, 'Biomass_', fileName, '.txt')) 
    end
end

