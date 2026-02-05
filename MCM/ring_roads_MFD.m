% clc;
addpath('./mcm');
clear variables;
% %close all;
seed = ceil(rand()*100);

numLayers = 1;
timeSteps = 10000;
% max number of packets that can jump from each node in 1 time step
tau = 1;

% load graph (edgelist format)
% fid = fopen("/Users/robertbenassai/Documents/UOC/k_shortest_" + ...
%     "path_betweenness/KSP/grapcd mcmhs_qspbwss/ring_roads/vx1/graph_variations/") ;
% data = textscan(fid,'%f %f %f','HeaderLines',1, 'delimiter',',');
% data = cell2mat(data);
% fclose(fid);
% 
% data(:, 1:2) = data(:, 1:2) + 1;
% read_A = spconvert(data);
% 
% [N1, N2] = size(read_A);
% W = sparse(zeros(max(N1, N2), max(N1, N2)));
% W(1:N1, 1:N2) = read_A;
% A = W+W';
% numNodes = length(A);
% %Structures that define the routing policy in the Shortest path routing
% [SPModelData.spBW,...
%         SPModelData.numExternalPacketRoutings,...
%         SPModelData.numLocalPacketRoutings,...
%         SPModelData.numStartsAtNode,...
%         SPModelData.numEndsAtNode]=...
%         SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap(A,numNodes,1);   
% SPModelData.spBW = sum(SPModelData.spBW,2);        
% %compute critical injection rate theorical
% %only accurate if the genration rate is the same for all node pairs
% rho_cs = tau*((numNodes-1)./(SPModelData.spBW + 2*(numNodes-1)));
% rho_c = min(rho_cs);   

% load the successor matrix directory. We'll have 1 successor matrix for 
% each node, which we'll store in 'succ'. E.g. if loaded succ_mat_i: start
% from column i+1 (due to matlab starting with 1) and each row indicates
% the final destination node. Thus, change columns (with row constant) as
% indicated in each cell to reach the final destination.

%tol_ls = [0, 30];
tol_ls = [0, 10, 20];
folders = dir('mcm/star_road_p_structure_vx2/');  % List all items in the 'mcm' directory
folders = folders([folders.isdir]);  % Keep only directories

% Remove '.' and '..' (current and parent directories)
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Filter directories containing 'ring'
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = length(ring_folders);

% delete folders with no path_ls file
% 
% ring_folders = folders(contains({folders.name}, 'ring_road_'));
% 
% keep = true(size(ring_folders));  % mask of which folders to keep
% 
% for k = 1:numel(ring_folders)
%     folderPath = fullfile(ring_folders(k).folder, ring_folders(k).name);
% 
%     % Look for files starting with "path_ls_mat"
%     files = dir(fullfile(folderPath, 'path_ls_mat*'));
% 
%     % If no such files → delete folder and mark for removal
%     if isempty(files)
%         fprintf('Deleting folder: %s\n', folderPath);
%         try
%             rmdir(folderPath, 's');   % delete recursively
%         catch ME
%             warning('Could not delete folder %s: %s', folderPath, ME.message);
%         end
%         keep(k) = false;
%     end
% end

% Update the folder list to only those that remain
% ring_folders = ring_folders(keep);

possGenRates = logspace(-2,log10(3),100);

nstat = zeros(length(possGenRates));
nstat_std = zeros(length(possGenRates));
noutstat = zeros(length(possGenRates));
noutstat_std = zeros(length(possGenRates));

numNodes = 100;
for tol_i = 1:length(tol_ls)
    DelivExp = zeros(length(possGenRates));
    DelivExp2 = zeros(length(possGenRates));
    DeltaNExp = zeros(numNodes, length(possGenRates)); % expected value of \Delta Q over all graphs
    DeltaNExp2 = zeros(length(possGenRates)); % expected value of \Delta Q**2
    % over all graphs

    tol = tol_ls(tol_i);
    h = waitbar(0, 'Initializing waitbar...');  
    wrong_files = 0;
    for file_i = 1:len_fold

        waitbar_label = sprintf("Processing several ring road graph" +  ...
            ' instances for tol %d. File %s', tol, ring_folders(file_i).name);
        waitbar(file_i / len_fold, ...
                h, waitbar_label)
        % str_mat = sprintf(['mcm/star_road_p_structure_vx2/%s/path_ls_' ...
        %     'mat_%d.mat'],ring_folders(file_i).name, max(tol_ls));

        str_mat = sprintf(['mcm/star_road_p_structure_vx1/%s/pen_path_ls_t_02_k_50_p_01' ...
            'mat_%d.mat'],ring_folders(file_i).name, max(tol_ls));
        
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

        for rhosIndex = 1:length(possGenRates)
            
            % setting node generation rates
            genRatePerMinuteGates = zeros(numNodes,1);
            genRatePerMinuteGates(:) = possGenRates(rhosIndex);       
        
            %define the "OD" matrix somehow, this is all to all but with the
            %given generationrate
            
            [SPModelData.possDest,SPModelData.possDestCumProb,SPModelData.possDestProb] =...
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
                    SPModelData.possDest,...
                    SPModelData.possDestCumProb,...
                    numNodes,1,timeSteps,seed); 
        
            % Add stationary number of packets and out rate
            tot_n = sum(numPackets, 1, "default");
            nstat(rhosIndex) = nstat(rhosIndex)+mean(tot_n(end-999:end));
            nstat_std(rhosIndex) = nstat_std(rhosIndex)+std(tot_n(end-999:end));
    
            noutstat(rhosIndex) = noutstat(rhosIndex)+mean(numPacketOutUnitTime(end-999:end));
            noutstat_std(rhosIndex) = noutstat_std(rhosIndex)+std(numPacketOutUnitTime(end-999:end));
            % compute microscopic vars obtained by the montecarlo
        
            normFactor = 1./sum(PExp);
            normFactor(normFactor == Inf) = 0;
            % Order parameter
            PExp = PExp*diag(normFactor);    
            pExp = packetsThatContinue./extractedPackets;
            DExp = mean(DExp,2);
            % compute queue increments obtained by the montecarlo
            sum_deltaNexp = 0;
            for j=1:numNodes
                p = polyfit(1:size(numPackets(j,:),2),numPackets(j,:),1);
                DeltaNExp(j,rhosIndex) = DeltaNExp(j,rhosIndex) + (p(1)); 
                sum_deltaNexp = sum_deltaNexp + (p(1));
            end    
                DeltaNExp2(rhosIndex) = DeltaNExp2(rhosIndex) + sum_deltaNexp*sum_deltaNexp;

            DelivExp(rhosIndex) = DelivExp(rhosIndex) + sum(numPacketOutUnitTime);%/sum(s_i, "all");
            DelivExp2(rhosIndex) = DelivExp2(rhosIndex) + sum(numPacketOutUnitTime)^2;
            packetsGenPerMinute(rhosIndex) = sum(genRatePerMinute);
            
        end
    end 

    % compute the average
    DelivExp = DelivExp/(len_fold-wrong_files);
    DeltaNExp = DeltaNExp/(len_fold-wrong_files);
    DelivExp2 = DelivExp2/(len_fold-wrong_files);
    DeltaNExp2 = DeltaNExp2/(len_fold-wrong_files);
    %computing stds
    DelivExp_std = sqrt(DelivExp2 - DelivExp.^2);
    DeltaNExp_std = sqrt(DeltaNExp2 - sum(DeltaNExp).^2);

    nstat = nstat/(len_fold-wrong_files);
    nstat_std = nstat_std/(len_fold-wrong_files);

    noutstat = noutstat/len_fold;
    noutstat_std = noutstat/len_fold;

    fname = sprintf('./sim_results/averaged/star/vx2/tol%d_mean_vx2.mat', tol);
    save(fname, 'DeltaNExp','DelivExp', 'DelivExp_std', ...
        'DeltaNExp_std', 'nstat','nstat_std','noutstat', ...
        'noutstat_std')
    close(h);
end

