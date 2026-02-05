% clc;
addpath('./mcm');
clear variables;
% %close all;
seed = ceil(rand()*100);

if isempty(gcp('nocreate'))
    parpool(8);
end

numLayers = 1;
timeSteps = 10000;
% max number of packets thacet can jump from each node in 1 time step
tau = 1;


tol_ls = [0, 10, 20];
% folders = dir('mcm/ring_road_p_structure_vx1/');  % List all items in the 'mcm' directory
% folders = dir('mcm/ring_road_p_structure_vx2/');  % List all items in the 'mcm' directory
folders = dir('mcm/star_road_p_structure_vx2/');  % List all items in the 'mcm' directory

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

% numNodes = 58;
numNodes = 100;

for tol_i = 1:length(tol_ls)
    DelivExp   = zeros(numRhos,1);
    DelivExp2  = zeros(numRhos,1);
    DeltaNExp  = zeros(numNodes,numRhos);
    DeltaNExp2 = zeros(numNodes,numRhos);   % if you really need it

    tol = tol_ls(tol_i);
    wrong_files = 0;
    for file_i = 1:len_fold

        file_i
        str_mat = sprintf('mcm/star_road_p_structure_vx2/%s/path_ls_mat_%d.mat', ...
            ring_folders(file_i).name, max(tol_ls));
        % str_mat = sprintf('mcm/ring_road_p_structure_vx2/%s/path_ls_mat_%d.mat', ...
        %     ring_folders(file_i).name, max(tol_ls));

        % 
        % str_mat = sprintf('mcm/ring_road_p_structure_vx1/%s/path_ls_mat_20.mat', ...
        %     ring_folders(file_i).name);
        
        %read cell array
        mat_struct = load(str_mat);
        fieldNames = fieldnames(mat_struct);

        if isempty(fieldNames)
            sprintf("entered %s\n", ring_folders(file_i).name);
            wrong_files = wrong_files+1;
            continue
        end
        PathLisArr = mat_struct.(fieldNames{1});

        if size(PathLisArr) == 0
            wrong_files = wrong_files+1;
            continue
        end

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

        parfor rhosIndex = 1:numRhos
            
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
        
            % monte carlo simulations adapted for tolerance bwss
            [numPackets,numPacketOutUnitTime,B_i_ext,e_i_ext,B_i,e_i,s_i,PExp,sigma,DExp,...
            extractedPackets,packetsThatContinue]...
                = cSPCongestion_statMem_dir_weighted_local_list_paths(...
                    genRatePerMinute,...
                    processingRatePerMinute,...
                    tol,...%SPModelData.DistRG,...
                    [],...%SPModelData.NumPathsRG,...
                    PathLisArr,...%path list cell array
                    possDest,...
                    possDestCumProb,...
                    numNodes,1,timeSteps,seed); 
              
            % --- statistics for this rho ---

            % Add stationary number of packets and out rate
            tot_n = sum(numPackets, 1, "default");

            nstat(rhosIndex) = mean(tot_n(end-999:end));
            nstat_std(rhosIndex) = std(tot_n(end-999:end));
    
            noutstat(rhosIndex) = mean(numPacketOutUnitTime(end-999:end));
            noutstat_std(rhosIndex) = std(numPacketOutUnitTime(end-999:end));
            % compute microscopic vars obtained by the montecarlo
        
            normFactor = 1./sum(PExp);
            normFactor(normFactor == Inf) = 0;
            % Order parameter
            PExp = PExp*diag(normFactor);    
            pExp = packetsThatContinue./extractedPackets;
            DExp = mean(DExp,2);

            % compute queue increments obtained by the montecarlo
            for j=1:numNodes
                p = polyfit(1:size(numPackets(j,:),2),numPackets(j,:),1);
                DeltaNExp(j,rhosIndex) = (p(1)); 
            end      

            DelivExp(rhosIndex) = sum(numPacketOutUnitTime);%/sum(s_i, "all");
            DelivExp2(rhosIndex) = sum(numPacketOutUnitTime)^2;

            packetsGenPerMinute(rhosIndex) = sum(genRatePerMinute);
        end
        DelivExp_std = sqrt(DelivExp2 - DelivExp.^2);
        % Write file for each graph

        s = string(ring_folders(file_i).name);
        nums = double(extract(s, digitsPattern));  % -> [123 2]
        graph_num = nums(1);

        % fname = sprintf('./sim_results/all_graphs/ring/vx1/tol%d_%d_vx1.mat', tol, graph_num);
        % fname = sprintf('./sim_results/all_graphs/ring/vx2/tol%d_%d_limq_vx2.mat', tol, graph_num);
        fname = sprintf('./sim_results/all_graphs/star/tol%d_%d_vx2.mat', tol, graph_num);

        save(fname, 'DeltaNExp','DelivExp', 'DelivExp_std', ...
        'nstat','nstat_std','noutstat', ...
        'noutstat_std')
        
    end 
end
