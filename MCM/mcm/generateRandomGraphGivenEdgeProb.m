%generateRandomGraphGivenEdgeProb
%   [ layerNetworkAdj ] = generateRandomGraphGivenEdgeProb( numNodes, connectivityThreshold, allowNotConnected )
%       numNodes :  number of nodes the network has
%       connectivityThreshold : probability of having an edge
%       layerNetworkAdj : adjacency matrix of the network, returned network
%           is connected. If no connected network can be found, "empty" is
%           returned
%       allowNotConnected : to remove connectedness restriction set to false.
%Albert Solé-Ribalta

function [ layerNetworkAdj,numTrials ] = generateRandomGraphGivenEdgeProb( numNodes, connectivityThreshold, allowNotConnected )

    if ~exist('allowNotConnected','var')
        allowNotConnected = false;
    end

    maxTrials = 1000;
    
    dataLayerNetwork = rand(1,(numNodes*(numNodes-1))/2)<=connectivityThreshold;
    layerNetwork = triu(ones(numNodes),1);
    layerNetwork(~~layerNetwork)=dataLayerNetwork;
    layerNetworkAdj = layerNetwork+layerNetwork';
    
    numTrials = 0;
    if(~allowNotConnected)
        [~, sizeConnComp]  = conncomp(graph(layerNetworkAdj));
        [~, numConnComp] = size(sizeConnComp);
        numTrials = 0;
        while ((numConnComp ~= 1) && (numTrials < maxTrials))         
            dataLayerNetwork = rand(1,(numNodes*(numNodes-1))/2)<=connectivityThreshold;
            layerNetwork = triu(ones(numNodes),1);
            layerNetwork(~~layerNetwork)=dataLayerNetwork;
            layerNetworkAdj = layerNetwork+layerNetwork';
            [~, sizeConnComp]  = conncomp(graph(layerNetworkAdj));
            [~, numConnComp] = size(sizeConnComp);

            numTrials = numTrials + 1;
        end

        if(numConnComp ~= 1)
            layerNetworkAdj = [];
        end
    end
end


% function [ layerNetworkAdj ] = generateRandomGraphGivenEdgeProb( numNodes, connectivityThreshold )
% 
%     maxTrials = 1000;
% 
%     dataLayerNetwork = rand(1,(numNodes*(numNodes-1))/2)<=connectivityThreshold;
%     layerNetwork = triu(ones(numNodes),1);
%     layerNetwork(~~layerNetwork)=dataLayerNetwork;
%     layerNetworkAdj = layerNetwork+layerNetwork';
%     layerNetworkLapl = diag(sum(layerNetwork+layerNetwork')) - (layerNetwork+layerNetwork');
%     secondEigVal = computeEigValVect(layerNetworkLapl,2);
%     
%     numTrials = 0;
%     while ((abs(secondEigVal)<1e-4) && (numTrials < maxTrials))         
%         dataLayerNetwork = rand(1,(numNodes*(numNodes-1))/2)<=connectivityThreshold;
%         layerNetwork = triu(ones(numNodes),1);
%         layerNetwork(~~layerNetwork)=dataLayerNetwork;
%         layerNetworkAdj = layerNetwork+layerNetwork';
%         layerNetworkLapl = diag(sum(layerNetwork+layerNetwork')) - (layerNetwork+layerNetwork');
%         secondEigVal = computeEigValVect(layerNetworkLapl,2);
%         numTrials = numTrials + 1;
%     end
% 
%     if(abs(secondEigVal)<1e-4)
%         layerNetworkAdj = [];
%     end
% end
% 
