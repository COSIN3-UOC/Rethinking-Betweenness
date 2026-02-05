
function [possDest,cumProbOfDest,probOfDest] = matrixOfAccesibleDestinations_congestionZone(gatesIndex, nodesIndex, genRatePerMinuteGates, numNodes)

    %each packet each node will receive will be returned to the gate
    for i=1:numNodes       
        possDest{i} = setdiff(1:numNodes,i);        
        probOfDest{i} = ones(numNodes-1,1) ./ (numNodes-1);
        cumProbOfDest{i} = cumsum(probOfDest{i})';
    end           
end
