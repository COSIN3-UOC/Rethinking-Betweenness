function [ DeltaN ] = computeDeltaN( A, P, p, D, I, numNodes )
    
    DeltaN = zeros(numNodes,1);

    for i=1:numNodes        
        DeltaN(i) = I(i) - D(i);
        for j=1:numNodes
            %if (i~=j)
            if(D(j)*P(i,j)*p(j) > 0)
                DeltaN(i) = DeltaN(i) + D(j)*P(i,j)*p(j);
            end
            %end
        end 
    end

end