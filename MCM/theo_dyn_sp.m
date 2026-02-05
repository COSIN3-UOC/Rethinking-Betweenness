% clc;
addpath('./mcm');
clear variables;
% %close all;
seed = ceil(rand()*100);

numLayers = 1;
timeSteps = 10000;
% max number of packets thacet can jump from each node in 1 time step
tau = 1;


tol_ls = [0];
% tol_ls = [0, 10, 20, 30];
% folders = dir('mcm/ring_road_p_structure_vx1_test/');  % List all items in the 'mcm' directory
folders = dir('mcm/ring_road_p_structure_vx2/');  % List all items in the 'mcm' directory

folders = folders([folders.isdir]);  % Keep only directories

% Remove '.' and '..' (current and parent directories)
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Filter directories containing 'ring'
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = length(ring_folders);
    
possGenRates = logspace(-2,log10(3),100);

numRhos = numel(possGenRates);

nstat        = zeros(numRhos,1);
nstat_std    = zeros(numRhos,1);
noutstat     = zeros(numRhos,1);
noutstat_std = zeros(numRhos,1);

numNodes = 58;

for tol_i = 1:length(tol_ls)

    tol = tol_ls(tol_i);
    h = waitbar(0, 'Initializing waitbar...');  
    wrong_files = 0;
    for file_i = 1:len_fold

        waitbar_label = sprintf("Processing several graph" +  ...
            " instances for tol %d. File %s", tol, ring_folders(file_i).name);
        waitbar(file_i / len_fold, ...
                h, waitbar_label)
        
        % load adjacency matrix
        g_num = double(extract(string(ring_folders(file_i).name), digitsPattern));
        fname_graph = sprintf(("/Users/robertbenassai/Documents/UOC/" + ...
            "k_shortest_path_betweenness/KSP/graphs_qspbwss/" + ...
            "ring_roads/vx2/graph_variations/ring%d.csv"), g_num);

        fid = fopen(fname_graph) ;
        data = textscan(fid,'%f %f %f','HeaderLines',1, 'delimiter',',');
        data = cell2mat(data);
        fclose(fid);
        
        data(:, 1:2) = data(:, 1:2) + 1;
        read_A = spconvert(data);
        
        [N1, N2] = size(read_A);
        W = sparse(zeros(max(N1, N2), max(N1, N2)));
        W(1:N1, 1:N2) = read_A;
        A = W+W';
        numNodes = length(A);

        % SPModelData contains:
        
        % - spBW: Shortest path betweenness
        % - numExternalPacketRoutings: IDK
        % - numLocalPacketRoutings: IDK
        % - numStartsAtNode: number of possible paths starting from each node:
        % nnodes-1
        % - numEndsAtNode: number of possible paths ending at each node: nnodes-1
        % - possDest: possible destinations (all nodes - starting node) for each
        % node
        % - possDestCumProb: cumulative probability of ending at each node (uniform)
        % - possDestProb: probability of ending at each node (uniform)
        
        % only to compute theoretical prediction:
        
        % - PredRG: the successor matrix
        % - PredProbRG: probabilities of each successor

        %possGenRates = linspace(0.01,3,100);

        congestedNodes = [];

        %Structures that define the routing policy
        [SPModelData.spBW,...
                SPModelData.numExternalPacketRoutings,...
                SPModelData.numLocalPacketRoutings,...
                SPModelData.numStartsAtNode,...
                SPModelData.numEndsAtNode]=...
                SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap(A,numNodes,1);
        
        SPModelData.spBW = sum(SPModelData.spBW,2); 
        
        %compute critical injection rate theorical
        %only accurate if the genration rate is the same for all node pairs
        rho_cs = tau*((numNodes-1)./(SPModelData.spBW + 2*(numNodes-1)));
        rho_c = min(rho_cs);   

        sp_bw = SPModelData.spBW;
        %%[SPModelData.DistRG,SPModelData.NumPathsRG,SPModelData.PredRG,SPModelData.PredProbRG,...
        [~,~,~,~,...
             SPModelData.srcMultRG,SPModelData.dstMultRG,SPModelData.probMultRG,...           
             ] = SPDataPathDegeneration_staticMemory(A,numNodes,1);  

        for rhosIndex = 1:numRhos
            
            % setting node generation rates
            genRatePerMinuteGates = zeros(numNodes,1);
            genRatePerMinuteGates(:) = possGenRates(rhosIndex);       
 
            %define the "OD" matrix somehow, this is all to all but with the
            %given generationrate
            
            [possDest,possDestCumProb,possDestProb] =...
                matrixOfAccesibleDestinations_AllToAll(...
                    find(genRatePerMinuteGates>0),...
                    find(genRatePerMinuteGates==0),...
                    genRatePerMinuteGates,...
                    numNodes);   
            genRatePerMinute = genRatePerMinuteGates;        
            processingRatePerMinute = ones(numNodes,1);    

            SPModelData.possDest = possDest;
            SPModelData.possDestProb = possDestProb;

            [DeltaNTeor(:,rhosIndex),rateNumPacketsDelivered,...
            B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs,...
            DTeor, pTeor, PTeor,congestedNodes] = SP_computeEtaGivenRho_dir_weigh_nonHomStart(...
            A,numNodes,numLayers,...
            genRatePerMinute,processingRatePerMinute,...
            SPModelData.possDest,SPModelData.possDestProb,SPModelData,congestedNodes);  
                          

        end
        % Write file for each graph

        s = string(ring_folders(file_i).name);
        nums = double(extract(s, digitsPattern));  % -> [123 2]
        graph_num = nums(1);

        fname = sprintf('./sim_results/all_graphs/penal/ring/vx2/tol%d_%d_vx2_theo.mat', tol, graph_num);
        % fname = sprintf('./sim_results/all_graphs/ring/vx2/tol%d_%d_vx2.mat', tol, graph_num);
        save(fname, 'DeltaNTeor', 'rho_cs', 'sp_bw')
        
    end 

    close(h);
end
