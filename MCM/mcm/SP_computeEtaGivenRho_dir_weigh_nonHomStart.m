function [DeltaN,rateNumPacketsDelivered, ...
    B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs,...
    D, p, P,congestedNodes] = ...
        SP_computeEtaGivenRho_dir_weigh_nonHomStart(A,numNodes,numLayers,...
        genRate,processingRate,possDest,possDestProb,...
        SPModelData, initialCongested)

    %%% INITIALIZATIONS %%%
    congestedNodes = [];
    if(~exist('initialCongested'))
        congestedNodes = [];
    else
        congestedNodes = initialCongested;
    end
    
    % Compute initial values of the model
    normFactorLocPackets = 1./sum(SPModelData.numLocalPacketRoutings);
    normFactorLocPackets(normFactorLocPackets == Inf) = 0;
    probTransitLocalPack = SPModelData.numLocalPacketRoutings*diag(normFactorLocPackets);    
    normFactorExtPackets = 1./sum(SPModelData.numExternalPacketRoutings);
    normFactorExtPackets(normFactorExtPackets == Inf) = 0;
    probTransitExternPack = SPModelData.numExternalPacketRoutings*diag(normFactorExtPackets);   
                    
    e_i_teor = SPModelData.numEndsAtNode;
    B_i_teor = SPModelData.spBW;
    
    D = ones(numNodes*numLayers,1);
    D(:) = min(genRate,processingRate);        
    D(congestedNodes) = processingRate(congestedNodes);
    
    p = B_i_teor./(B_i_teor + e_i_teor);   
    p(p~=p) = 0;
    assert(all(p~=Inf));       
    
    P = probTransitLocalPack + probTransitExternPack;     
    normFactor = 1./sum(P);
    normFactor(normFactor == Inf) = 0;
    P = P*diag(normFactor);    
    % End compute initial values of the model
        
    %%% RUN THE MODEL %%%    
    congestedNodes = zeros(numNodes,1);
    congestedNodes(initialCongested) = 1;
    newCongestedDetected = true;
    while(newCongestedDetected)
        notCongestedIndices = find(congestedNodes==0);
        congestedNodesIndices = find(congestedNodes==1);
        
        [ D, p, P, B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs] = ...
            sp_computeMCMParams_dir_weighted_dynObsBwComp_C( ...
            congestedNodesIndices, A, D, P, p , ...
            genRate, processingRate, possDest, possDestProb, ...
            numNodes,numLayers);            
        
        [maxDValue,maxDNode] = max(D(notCongestedIndices));

        if(maxDValue > processingRate(maxDNode))
            congestedNodes(notCongestedIndices(maxDNode)) = 1;
            newCongestedDetected = true;
        else
            newCongestedDetected = false;
        end
    end         

    DeltaN = computeDeltaN( A, P, p, D, genRate, numNodes*numLayers );
%    rateNumPacketsDelivered(rhoIndex) = sum(DTeor.*(1-pTeor));
    rateNumPacketsDelivered = [];
    congestedNodes = find(congestedNodes == 1);
end

function [ D, p, P, B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs] = sp_computeMCMParams_dir_weighted_dynObsBwComp_C( ...
    congestedNodes, A, D, P, p, genRate, processingRate, possDest, possDestProb,...
    numNodes, numLayers)
    
    T = 1000;
        
    congestedVector = zeros(numNodes*numLayers,1);
    congestedVector(congestedNodes) = 1;
            
    [P,D,p, B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs,~,~,~] = ...
        sp_computeMCMParams_dir_weighted_dynObsBwComp_C_aux( A, D, P, p, congestedVector,...
        genRate, processingRate, possDest,possDestProb, numNodes,numLayers,T);    
end
