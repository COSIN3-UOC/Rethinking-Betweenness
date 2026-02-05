// to compile use: mex -largeArrayDims sp_computeMCMParams_dir_weighted_dynObsBwComp_C_aux.cpp -I/usr/local/Cellar/boost/1.80.0/include/
#define TRUE 1
#define FALSE 0
#define DEBUG 0
#include "mex.h"
#include <string.h>
#include <math.h>
#define MAX_COLA 2000
#include "predQueueCircStatic.h"
#include <boost/heap/fibonacci_heap.hpp>
#define MAX_VAL 1000000.0
#define SAT_DATA_SIZE 2000

char outputString[100];

//to use the heaps
struct dk_Node{
    int nodeId;
    double distance;
    dk_Node(int id, double dist) : nodeId(id),distance(dist) { } //<-- this is a constructor
};
struct compare_dk_Node{
    bool operator()(const dk_Node& n1, const dk_Node& n2) const{
        return n1.distance > n2.distance;
    }
};

void SPEdgeNodeBetweennessWeightNetCong(
        const mxArray * A, int numNodes, int numLayers, double * inOutFluxRatios, 
        const mxArray * possDest, const mxArray * possDestProb,double * genRate,
        double * spNBW, double * spEBW_ExcludeInit, double * spEBW_Init, double * numEndAtNode, double * numStartAtNode);
void printSquareMatrix(double * P, int sizeSupra);
void printMatrix(double * P, int rows, int cols);
void printMatrixInt(int * P, int rows, int cols);
void projectRowsWithAddition(double * resultingVect, double * matrix, int numRows, int numColumns);
void projectColumnsWithAddition(double * resultingVect, double * matrix, int numRows, int numColumns);
void printMatrixWithIndex(double * P, int rows, int cols);
void fillDesinationsProbMatrix(
        double * possDestProbMatrix,
        const mxArray * possDest,
        const mxArray * possDestProb, int sizeSupra);
void obtain_PNormConst(
        double cumProbTransitLocalPack,
        double cumProbTransitExternPack,
        double * normConst,
        double *oneMinusNormConst);
void recomputeJumpingPorbabilities(
        double * B_i_obs,
        double * e_i_obs,
        double * EB_Ext_obs,
        double * EB_Int_obs,
        double * pCondExtPacket, 
        double * probTransitExternPack,
        double * probTransitLocalPack,
        int sizeSupra);

//gloval vars for the model
double M_sigma[SAT_DATA_SIZE]; 
double M_inOutFluxRatios[SAT_DATA_SIZE];
double M_s_i_obs[SAT_DATA_SIZE];
double M_totalNumRoutings[SAT_DATA_SIZE];
//global vars for the betweenness
staticQueue BW_S;
staticQueue BW_SForward;   
double BW_visited[SAT_DATA_SIZE];
double BW_sigma[SAT_DATA_SIZE];
double BW_d[SAT_DATA_SIZE]; 
staticQueue BW_P[SAT_DATA_SIZE];
staticQueue BW_Successors[SAT_DATA_SIZE];   
double BW_delta[SAT_DATA_SIZE];             
double BW_possDestProbMatrix[SAT_DATA_SIZE*SAT_DATA_SIZE];           
double cumProbTransitLocalPack[SAT_DATA_SIZE];
double cumProbTransitExternPack[SAT_DATA_SIZE];
double BW_p[SAT_DATA_SIZE];

double spNBW_Temp[SAT_DATA_SIZE];   
double spEBW_ExcludeInit_temp[SAT_DATA_SIZE*SAT_DATA_SIZE];   
double spEBW_Init_temp[SAT_DATA_SIZE*SAT_DATA_SIZE];     
double numEndAtNode_temp[SAT_DATA_SIZE];      
double numStartAtNode_temp[SAT_DATA_SIZE];     

double normalizedD[SAT_DATA_SIZE];
double normalizedP[SAT_DATA_SIZE*SAT_DATA_SIZE];
double normalizedp[SAT_DATA_SIZE];
double normalizedDOld[SAT_DATA_SIZE];
double normalizedPOld[SAT_DATA_SIZE*SAT_DATA_SIZE];
double normalizedpOld[SAT_DATA_SIZE];

//to store the nodes and the distance to them
boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> > dk_heap;        
boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> >::handle_type heapHandles[MAX_COLA];    

    
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){       
    
    //vars
    int t;
    int i;
    int j;
    int sizeSupra;
    int converged = FALSE;
    double normD = 0;
    double normP = 0;
    double normp = 0;   
    double threshold = 1e-6;
    double normConst = 0;
    double oneMinusNormConst = 0;
    double anyLowerThanOne;
    
    //check inputs
    if (mxGetM( prhs[0] ) != mxGetN( prhs[0] )) mexErrMsgTxt( "Input matrix G needs to be square." );
    if(mxIsSparse(prhs[0])==0) mexErrMsgTxt( "Distance Matrix must be sparse" );    
    
    //input data
    double * DOld = mxGetPr(prhs[1]);
    double * POld = mxGetPr(prhs[2]);
    double * pOld = mxGetPr(prhs[3]);
    double * congestedNodes = mxGetPr(prhs[4]);
    //double * numEndAtNode = mxGetPr(prhs[5]);
    //double * BW_Ext = mxGetPr(prhs[6]);    
//    double * pCondExtPacket = mxGetPr(prhs[5]);  //(double *)mxMalloc(sizeof(double)*sizeSupra);
//    double * probTransitLocalPack = mxGetPr(prhs[6]);  //(double *)mxMalloc(sizeof(double)*sizeSupra*sizeSupra);//mxGetPr(prhs[6]);  
//    double * probTransitExternPack = mxGetPr(prhs[7]);  //(double *)mxMalloc(sizeof(double)*sizeSupra*sizeSupra);//mxGetPr(prhs[7]);      
//     double * pCondExtPacket = (double *)mxMalloc(sizeof(double)*sizeSupra);
//     double * probTransitLocalPack = (double *)mxMalloc(sizeof(double)*sizeSupra*sizeSupra);//mxGetPr(prhs[6]);  
//     double * probTransitExternPack = (double *)mxMalloc(sizeof(double)*sizeSupra*sizeSupra);//mxGetPr(prhs[7]);  
    
    double * genRate = mxGetPr(prhs[5]);    
    double * processingRate = mxGetPr(prhs[6]);    
    const mxArray * possDest = prhs[7];
    const mxArray * possDestProb = prhs[8];            
    int numNodes = (int)mxGetPr(prhs[9])[0];
    int numLayers = (int)mxGetPr(prhs[10])[0];    
    int T = (int)mxGetPr(prhs[11])[0]; //maximum time steps until convergence
    sizeSupra = numNodes*numLayers;        
    
    //sprintf(outputString,"fprintf('External\\n');");
    //mexEvalString(outputString);  
    //printMatrixWithIndex(&cumProbTransitExternPack[0], 1, sizeSupra);
    //sprintf(outputString,"fprintf('Local\\n');",t,normD,normP,normp);
    //mexEvalString(outputString);  
    //printMatrixWithIndex(&cumProbTransitLocalPack[0], 1, sizeSupra);    
    //return;
    
    //output data
    plhs[0] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );
    double * P = mxGetPr(plhs[0]);    
    memset ( P, 0, sizeof(double)*sizeSupra*sizeSupra);     
    plhs[1] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );
    double * D = mxGetPr(plhs[1]);    
    plhs[2] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );
    double * p = mxGetPr(plhs[2]);            
    plhs[3] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );        
    double * B_i_obs = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );    
    double * e_i_obs = mxGetPr(plhs[4]);    
    plhs[5] = mxCreateDoubleMatrix( sizeSupra,sizeSupra, mxREAL );    
    double * EB_Ext_obs = mxGetPr(plhs[5]);
    plhs[6] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );    
    double * EB_Int_obs = mxGetPr(plhs[6]);          
    plhs[7] = mxCreateDoubleMatrix( sizeSupra, 1, mxREAL );    
    double * pCondExtPacket = mxGetPr(plhs[7]);    
    plhs[8] = mxCreateDoubleMatrix( sizeSupra,sizeSupra, mxREAL );    
    double * probTransitLocalPack = mxGetPr(plhs[8]);
    plhs[9] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );    
    double * probTransitExternPack = mxGetPr(plhs[9]);        

    fillDesinationsProbMatrix(&BW_possDestProbMatrix[0], possDest, possDestProb, sizeSupra);    
    
    //this is to set the values in case there is no ratio lower than 1     
    for(i = 0 ; i < sizeSupra ; i++) M_inOutFluxRatios[i] = 1;                   
    SPEdgeNodeBetweennessWeightNetCong(
            prhs[0],numNodes,numLayers,
            &M_inOutFluxRatios[0],possDest,
            possDestProb, genRate,
            B_i_obs, EB_Ext_obs, EB_Int_obs, e_i_obs, M_s_i_obs);        
    
    recomputeJumpingPorbabilities(
            B_i_obs,
            e_i_obs,
            EB_Ext_obs,
            EB_Int_obs,
            pCondExtPacket,
            probTransitExternPack,
            probTransitLocalPack,
            sizeSupra);
    
    //compute the cummulative probability of jumping, we are interested if some of the columns sum zero
    projectColumnsWithAddition(&cumProbTransitLocalPack[0], probTransitLocalPack, sizeSupra, sizeSupra);
    projectColumnsWithAddition(&cumProbTransitExternPack[0], probTransitExternPack, sizeSupra, sizeSupra);

      //sprintf(outputString,"fprintf('EB_Int_obs\\n');"); mexEvalString(outputString);
      //printMatrix(&EB_Int_obs[0],sizeSupra,sizeSupra);
    //return;
//      sprintf(outputString,"fprintf('EB_Ext_obs\\n');"); mexEvalString(outputString);
//      printMatrix(&EB_Ext_obs[0],sizeSupra,sizeSupra);
//      sprintf(outputString,"fprintf('EB_Int_obs\\n');"); mexEvalString(outputString);
//      printMatrix(&EB_Int_obs[0],sizeSupra,sizeSupra);
//      sprintf(outputString,"fprintf('B_i_obs\\n');"); mexEvalString(outputString);
//      printMatrix(&B_i_obs[0],1,sizeSupra);
//      sprintf(outputString,"fprintf('e_i_obs\\n');"); mexEvalString(outputString);
//      printMatrix(&e_i_obs[0],1,sizeSupra);
//      sprintf(outputString,"fprintf('M_s_i_obs\\n');"); mexEvalString(outputString);
//      printMatrix(&M_s_i_obs[0],1,sizeSupra);
//      return;    
        
    memset( M_sigma, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( M_inOutFluxRatios, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( M_s_i_obs, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( M_totalNumRoutings, 0, sizeof(double)*SAT_DATA_SIZE); 

    memset( normalizedD, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( normalizedP, 0, sizeof(double)*SAT_DATA_SIZE*SAT_DATA_SIZE);
    memset( normalizedp, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( normalizedDOld, 0, sizeof(double)*SAT_DATA_SIZE); 
    memset( normalizedPOld, 0, sizeof(double)*SAT_DATA_SIZE*SAT_DATA_SIZE);
    memset( normalizedpOld, 0, sizeof(double)*SAT_DATA_SIZE); 
    
    //normalize the vectors p, P and D
    normD = 0;
    normP = 0;
    normp = 0;        
    for(i=0;i<sizeSupra;i++){
        normD = normD + DOld[i];
        normp = normp + pOld[i];
        for(j=0;j<sizeSupra;j++){                    
            normP = normP + POld[i + j*sizeSupra];
        }
    }        
    for(i=0;i<sizeSupra;i++){
        normalizedDOld[i] = D[i]/normD;
        normalizedpOld[i] = p[i]/normp;
        for(j=0;j<sizeSupra;j++){                    
            normalizedPOld[i + j*sizeSupra] = P[i + j*sizeSupra]/normP;
        }
    }        
    //End normalize the vectors p, P and D

    t = 0;
    while (!converged && t < T ){
         
        //compute the sigmas        
        memset(&M_sigma[0],0,sizeof(double)*sizeSupra);        
        for(i=0;i<sizeSupra;i++){
            for(j=0;j<sizeSupra;j++){
                if(i!=j){
                    //M_sigma[i] = M_sigma[i] + A[i+j*sizeSupra]*POld[i+j*sizeSupra]*DOld[j]*pOld[j]; 
                    M_sigma[i] = M_sigma[i] + POld[i+j*sizeSupra]*DOld[j]*pOld[j];  // multiplication for A not needed, same info 
                                                                                // encoded in P
                }
            }        
        }
        
//         sprintf(outputString,"fprintf('M_sigma\\n');"); mexEvalString(outputString);         
//         printMatrix(M_sigma,1,numNodes);
//         return;
//         
        memset(D,0,sizeof(double)*sizeSupra);        
        for(i=0;i<sizeSupra;i++){
            if(congestedNodes[i]==0){
                D[i] = genRate[i] + M_sigma[i];            
            }else{
                D[i] = processingRate[i];
            }     
        }

        ////////RECOMPUTE p and P//////////
        //compute in out flux ratios
        anyLowerThanOne = 0;
        for(i = 0 ; i < sizeSupra ; i++){
            M_inOutFluxRatios[i] = processingRate[i]/(M_sigma[i]+genRate[i]);
            if(M_inOutFluxRatios[i] < 1){
                anyLowerThanOne = 1;                
            }
            //sprintf(outputString,"fprintf('inOutFluxRatios[%i] = %f\\n');",i,M_inOutFluxRatios[i]); mexEvalString(outputString);                   
        }        
        
        if(anyLowerThanOne==1){                

            SPEdgeNodeBetweennessWeightNetCong(
                prhs[0], numNodes, numLayers, &M_inOutFluxRatios[0], possDest, possDestProb, genRate,
                B_i_obs, EB_Ext_obs, EB_Int_obs, e_i_obs, &M_s_i_obs[0]);                

            //sprintf(outputString,"fprintf('EB_Ext_obs\\n');");mexEvalString(outputString);
            //printMatrix(EB_Ext_obs,sizeSupra,sizeSupra);            
            //sprintf(outputString,"fprintf('EB_Int_obs\\n');");mexEvalString(outputString);
            //printMatrix(EB_Int_obs,sizeSupra,sizeSupra);            
                        
            recomputeJumpingPorbabilities(
                    B_i_obs,
                    e_i_obs,
                    EB_Ext_obs,
                    EB_Int_obs,
                    pCondExtPacket,
                    probTransitExternPack,
                    probTransitLocalPack,
                    sizeSupra);

        }        
        //sprintf(outputString,"fprintf('probTransitExternPack\\n');");mexEvalString(outputString);
        //printMatrix(probTransitExternPack,sizeSupra,sizeSupra);            
        //sprintf(outputString,"fprintf('probTransitLocalPack\\n');");mexEvalString(outputString);
        //printMatrix(probTransitExternPack,sizeSupra,sizeSupra);                    
        ////////END RECOMPUTE p //////////               

        
        memset(p,0,sizeof(double)*sizeSupra);        
        for(i=0;i<sizeSupra;i++){
            if((M_sigma[i]+genRate[i]) > 0){
                p[i] = pCondExtPacket[i]*(M_sigma[i]/(M_sigma[i]+genRate[i])) +  genRate[i]/(M_sigma[i]+genRate[i]);            
            }else{
                p[i] = 0;
            }
        }
                
        memset(P,0,sizeof(double)*sizeSupra*sizeSupra);
        for(i=0;i<sizeSupra;i++){
            normConst = (genRate[i]/(pCondExtPacket[i]*M_sigma[i]+genRate[i]));
            oneMinusNormConst = (pCondExtPacket[i]*M_sigma[i]/(pCondExtPacket[i]*M_sigma[i]+genRate[i]));
            
            //adapt the norm constants in case one of the columns is full of zeros
            //obtain_PNormConst(cumProbTransitLocalPack[i],cumProbTransitExternPack[i],&normConst,&oneMinusNormConst);                
                
            for(j=0;j<sizeSupra;j++){       
                if( (pCondExtPacket[i]*M_sigma[i]+genRate[i]) > 0 ){                                                                                                    
                    P[j + i*sizeSupra] = normConst*probTransitLocalPack[j + i*sizeSupra] 
                            + 
                        oneMinusNormConst*probTransitExternPack[j + i*sizeSupra];   
                }else{
                    P[j + i*sizeSupra] = 0;
                }
                //sprintf(outputString,"fprintf('PExp[%i,%i] = %f\\n');",j+1,i+1,P[j + i*sizeSupra]);mexEvalString(outputString);
                //sprintf(outputString,"fprintf('\\tprobTransitLocalPack[%i,%i] = %f\\n');",j+1,i+1,probTransitLocalPack[j + i*sizeSupra]);mexEvalString(outputString);
                //sprintf(outputString,"fprintf('\\tprobTransitExternPack[%i,%i] = %f\\n');",j+1,i+1,probTransitExternPack[j + i*sizeSupra]);mexEvalString(outputString);
            }        
        }

        //sprintf(outputString,"fprintf('PExp\\n');");mexEvalString(outputString);
        //printMatrix(P,sizeSupra,sizeSupra);            
        //return;
        
        //increment t
        t = t + 1;
        
        //normalize the vectors p, P and D
        normD = 0;
        normP = 0;
        normp = 0;        
        for(i=0;i<sizeSupra;i++){
            normD = normD + D[i];
            normp = normp + p[i];
            for(j=0;j<sizeSupra;j++){                    
                normP = normP + P[i + j*sizeSupra];
            }
        }        
        for(i=0;i<sizeSupra;i++){
            normalizedD[i] = D[i]/normD;
            normalizedp[i] = p[i]/normp;
            for(j=0;j<sizeSupra;j++){                    
                normalizedP[i + j*sizeSupra] = P[i + j*sizeSupra]/normP;
            }
        }        
                                
        //compute the norm of the matrix
        normD = 0;
        normP = 0;
        normp = 0;        
        for(i=0;i<sizeSupra;i++){
            normD = normD + pow(normalizedDOld[i] - normalizedD[i],2);
            normp = normp + pow(normalizedpOld[i] - normalizedp[i],2);
            for(j=0;j<sizeSupra;j++){                    
                normP = normP + pow(normalizedPOld[i + j*sizeSupra] - normalizedP[i + j*sizeSupra],2);
            }
        }        
        normD = sqrt(normD);
        normP = sqrt(normP);
        normp = sqrt(normp);        

        
        if((normD < threshold) & (normP < threshold) & (normp < threshold)){
        //if((normP < threshold) & (normp < threshold)){
            converged = TRUE;
        }else{
            memcpy(pOld,p,sizeof(double)*sizeSupra);
            memcpy(DOld,D,sizeof(double)*sizeSupra);
            memcpy(POld,P,sizeof(double)*sizeSupra*sizeSupra);
            
            memcpy(normalizedpOld,normalizedp,sizeof(double)*sizeSupra);
            memcpy(normalizedDOld,normalizedD,sizeof(double)*sizeSupra);
            memcpy(normalizedPOld,normalizedP,sizeof(double)*sizeSupra*sizeSupra);
            
        }

        
        
        //if(DEBUG){
            sprintf(outputString,"fprintf('it %i normD = %f, normP = %f, normp = %f\\n');",t,normD,normP,normp);
            mexEvalString(outputString);
        //}
    }
        
//    mxFree(pCondExtPacket);
//    mxFree(probTransitLocalPack);
//    mxFree(probTransitExternPack);   
        
//     mxFree(M_sigma);
//     mxFree(M_inOutFluxRatios);
//     mxFree(M_s_i_obs); 
//     mxFree(M_totalNumRoutings);
    
}     

void obtain_PNormConst(double cumProbTransitLocalPack, double cumProbTransitExternPack, double * normConst,double *oneMinusNormConst){
    double tolerance = 0.000000001;        
    
    if(cumProbTransitLocalPack <= tolerance && cumProbTransitExternPack <= tolerance){
        *normConst = 0;
        *oneMinusNormConst = 0;
    }
    if(cumProbTransitLocalPack <= tolerance && cumProbTransitExternPack > tolerance){
        *normConst = 0;
        *oneMinusNormConst = 1;
    }
    if(cumProbTransitLocalPack > tolerance && cumProbTransitExternPack <= tolerance){
        *normConst = 1;
        *oneMinusNormConst = 0;
    }            
}


// SPEdgeNodeBetweennessWeightNet(
//     prhs[0], numNodes, numLayers, inOutFluxRatios, possDest, possDestProb, genRate,
//     B_i_obs, EB_Ext_obs, EB_Int_obs, e_i_obs, s_i_obs);        

void SPEdgeNodeBetweennessWeightNetCong(
        const mxArray * A, int numNodes, int numLayers, double * inOutFluxRatios, 
        const mxArray * possDest, const mxArray * possDestProb,double * genRate,
        double * spNBW, double * spEBW_ExcludeInit, double * spEBW_Init, double * numEndAtNode, double * numStartAtNode){

    double tolerance = 0.000000001;    
    int sizeSupra = numNodes*numLayers;    
    double aux;
    double *Gpr = mxGetPr(A);
    size_t *Gir = mxGetIr(A);
    size_t *Gjc = mxGetJc(A);   
    //double * GData = mxGetData(prhs[0]);
    double edgeWeight = 0;
    double min_dist = -1;
    
    long int sIndEdges;
    long int eIndEdges;    
    long int nInd;
    int v,s,w,w2,l,temp,t,min_v;
      
    double numPacketReachDest;
                   
    //double * BW_visited = (double *)mxMalloc(sizeof(double)*sizeSupra);     
    //double * BW_sigma = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    //double * BW_sigmaOrig = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    //double * BW_d = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    //staticQueue * S = (staticQueue *)mxMalloc(sizeof(staticQueue));
    //staticQueue * P = (staticQueue *)mxMalloc(sizeof(staticQueue)*sizeSupra);
    //double * BW_delta = (double *)mxMalloc(sizeof(double)*sizeSupra);             
    //double * BW_possDestProbMatrix = (double *)mxMalloc(sizeof(double)*sizeSupra*sizeSupra);           
    //initialze the tables, this can maybe improved using a binaySearchTree
    //fillDesinationsProbMatrix(&BW_possDestProbMatrix[0], possDest, possDestProb, sizeSupra);
        
    //reinit all data, just in case!!
    memset ( spNBW, 0, sizeof(double)*sizeSupra ); 
    memset ( spEBW_ExcludeInit, 0, sizeof(double)*sizeSupra*sizeSupra ); 
    memset ( spEBW_Init, 0, sizeof(double)*sizeSupra*sizeSupra ); 
    memset ( numEndAtNode, 0, sizeof(double)*sizeSupra ); 
    memset ( numStartAtNode, 0, sizeof(double)*sizeSupra );                 
        
    //compute the percentage of path that will go out of the node
    for(w = 0; w<sizeSupra; w++){
        if(inOutFluxRatios[w]>=1)
            inOutFluxRatios[w] = 1;
    }    
    
    for( s = 0 ; s<numNodes ; s++){ 
    
        //sprintf(outputString,"fprintf('Source = %i, numNodes = %i, sizeSupra = %i, Inf = %f\\n');",s+1,numNodes,sizeSupra,MAX_VAL);
        //mexEvalString(outputString);

        //initializations
        memset ( spNBW_Temp, 0, sizeof(double)*sizeSupra ); 
        memset ( spEBW_ExcludeInit_temp, 0, sizeof(double)*sizeSupra*sizeSupra ); 
        memset ( spEBW_Init_temp, 0, sizeof(double)*sizeSupra*sizeSupra ); 
        memset ( numEndAtNode_temp, 0, sizeof(double)*sizeSupra ); 
        memset ( numStartAtNode_temp, 0, sizeof(double)*sizeSupra );                             
        create(&BW_S);        
        create(&BW_SForward); 
        for( w = 0 ; w<sizeSupra ; w++){            
            create(&BW_P[w]);     
            create(&BW_Successors[w]);             
            BW_sigma[w] = 0;                        
            BW_d[w] = MAX_VAL;
            BW_visited[w] = 0;
            
            //store the handle in a table then use the table to access the handle
            heapHandles[w] = dk_heap.push(dk_Node(w,BW_d[w]));                         
        }            
        BW_d[s] = 0;
        (*heapHandles[s]).distance = 0; dk_heap.decrease(heapHandles[s]);        
        BW_sigma[s] = 1;                
        //end initializations

        min_v = dk_heap.top().nodeId;
        min_dist = dk_heap.top().distance;
        dk_heap.pop();        
                        
        //while( min_dist < mxGetInf()){
        while( min_dist < MAX_VAL){
        //Min_V := S;
        //Min_Dist := 0.0;
        //while Min_Dist /= Plus_Infinity loop

            
            v = min_v;   
            BW_visited[v] = 1;
            add(&BW_S,v);                                                     
            //V := Min_V;
            //Visited(V) := True;
            //Push(V, St);
            
            if(DEBUG){
                sprintf(outputString,"fprintf('Exploring %i\\n');",v+1);
                mexEvalString(outputString);                                                     
            }            
            
            //for each neigbor of v
            sIndEdges = Gjc[v];
            eIndEdges = Gjc[v+1] - 1;            

            if(DEBUG){
                sprintf(outputString,"fprintf('min_v = %i, min_dist = %f\\n');",min_v+1, min_dist);
                mexEvalString(outputString);                                                                         
            }
                
            for( nInd = sIndEdges ; nInd<=eIndEdges ; nInd++ ){                        
            //while Has_Next(El) loop
                w = Gir[nInd];
                edgeWeight = Gpr[nInd];

                if(DEBUG){
                    sprintf(outputString,"fprintf('\tEdge from w=%i to w=%i and edgeWeight = %f\\n');",v+1, w+1, edgeWeight);
                    mexEvalString(outputString);                                                                         
                }                    
                
                if(DEBUG){                    
                    sprintf(outputString,"fprintf('\\tedges out %i\\n');",w);
                    mexEvalString(outputString);     
                }

                if (v != w){
                //if V /= W then                    
                    if(BW_d[w] > BW_d[v] + edgeWeight){ // new path to the vertex
                        BW_d[w] = BW_d[v] + edgeWeight;
                        (*heapHandles[w]).distance = BW_d[w]; dk_heap.decrease(heapHandles[w]);                                                
                        BW_sigma[w] = BW_sigma[v];
                        create(&BW_P[w]);
                        add(&BW_P[w],v);
                        
                        if(DEBUG){   
                            sprintf(outputString,"fprintf('\\t\\tnew path to vertex %i with dist  = %f\\n');",w+1, BW_d[w]);
                            mexEvalString(outputString);                                                     
                        }
                    }else if( fabs(BW_d[w] - BW_d[v] - edgeWeight) < tolerance){ // the current path is of equal length
                        BW_sigma[w] = BW_sigma[w] + BW_sigma[v];
                        add(&BW_P[w],v);                                                           
                    
                        if(DEBUG){                           
                            sprintf(outputString,"fprintf('\\t\\tpath to equal length to %i , path difference  =%f\\n');",w+1, fabs(BW_d[w] - BW_d[v] - edgeWeight));
                            mexEvalString(outputString);          
                            sprintf(outputString,"fprintf('\\t\\told length = %f, dist to %i = %f and edgeweight = %f\\n');",BW_d[w],v+1,BW_d[v],edgeWeight);
                            mexEvalString(outputString);                             
                        }
                    }                                                            
                 }
                 //end if;
            }
            //end loop;                       
            if(!dk_heap.empty()){
                min_v = dk_heap.top().nodeId;
                min_dist = dk_heap.top().distance;
                dk_heap.pop();                    
            }else{
                min_v = -1;
                min_dist = mxGetInf();            
            }            
        }                        
        
        //empty the heap if there is still some node
        while(!dk_heap.empty()){
            dk_heap.pop();                            
        }
        
        //keep an original copy of the BW_sigma
        //memcpy(BW_sigmaOrig,BW_sigma,sizeof(double)*sizeSupra);                
        for( w = 0 ; w<sizeSupra ; w++){
            BW_delta[w] = 0;            
        }          
        
                
        while(!empty(&BW_S)){
            w = pop(&BW_S);
            add(&BW_SForward,w);         
            if(DEBUG){
                sprintf(outputString,"fprintf('Distance from (s=%i) to (w=%i) = %f\\n');",s+1,w+1,BW_d[w]);
                mexEvalString(outputString);                                                         
            }
                
            while(!empty(&BW_P[w])){
                v = pop(&BW_P[w]);
                add(&BW_Successors[v],w);
                
                //compute the amount of paths that will start at w (the "1" of (1 + \delta_w))
                numPacketReachDest = genRate[s]*BW_possDestProbMatrix[w*sizeSupra + s];                
                //sprintf(outputString,"fprintf('\\tPredecessor of %i -> %i, numPackets = %f\\n');",v+1,w+1,numPacketReachDest);
                //mexEvalString(outputString);   
                BW_delta[v] = BW_delta[v] + (BW_sigma[v]/BW_sigma[w])*(numPacketReachDest + BW_delta[w]);
                if(v != s){
                    spEBW_ExcludeInit_temp[v*sizeSupra + w] = 
                            spEBW_ExcludeInit_temp[v*sizeSupra + w] + (BW_sigma[v]/BW_sigma[w])*(numPacketReachDest + BW_delta[w]);
                }                
                if(v == s){
                    spEBW_Init_temp[v*sizeSupra + w] = spEBW_Init_temp[v*sizeSupra + w] + (BW_sigma[v]/BW_sigma[w])*(numPacketReachDest + BW_delta[w]);
                    numStartAtNode_temp[v] = 
                           numStartAtNode_temp[v] + (BW_sigma[v]/BW_sigma[w])*(numPacketReachDest + BW_delta[w]);
                }                
                numEndAtNode_temp[w] = numEndAtNode_temp[w] + (BW_sigma[v]/BW_sigma[w])*numPacketReachDest;
            }
                
            if(w!=s){                                
                spNBW_Temp[w] = spNBW_Temp[w] + BW_delta[w];
            }
        }
        
        //sprintf(outputString,"fprintf('spNBW_Temp\\n');"); mexEvalString(outputString);  
        //printMatrix(spNBW_Temp,1,sizeSupra);
        //sprintf(outputString,"fprintf('genRate\\n');"); mexEvalString(outputString);  
        //printMatrix(genRate,1,sizeSupra);
        //sprintf(outputString,"fprintf('\\t\\tBW_sigma\\n');");mexEvalString(outputString);                                                     
        //printMatrix(BW_sigma,1,sizeSupra);
        //sprintf(outputString,"fprintf('numEndAtNode_temp\\n');"); mexEvalString(outputString);  
        //printMatrix(numEndAtNode_temp,1,sizeSupra);        
        
        // START: forward correction of the betweenness
        memset(BW_p,0,SAT_DATA_SIZE); 
        for( w = 0 ; w<sizeSupra ; w++){
            BW_delta[w] = 0;            
            if(w != s){
                if((spNBW_Temp[w] + numEndAtNode_temp[w]) > 0){
                    BW_p[w] = (spNBW_Temp[w]/(spNBW_Temp[w] + numEndAtNode_temp[w]));
                }else{
                    BW_p[w] = 0;
                }
            }else{
                BW_p[w] = 1;
            }
        }       
        BW_delta[s] = numStartAtNode_temp[s];
        numStartAtNode[s] = numStartAtNode[s] + inOutFluxRatios[s]*numStartAtNode_temp[s];

//         sprintf(outputString,"fprintf('BW_p\\n');"); mexEvalString(outputString);        
//         printMatrix(BW_p,1,sizeSupra);
//         return;
        
        while(!empty(&BW_SForward)){
            v = pop(&BW_SForward);            

            while(!empty(&BW_Successors[v])){
                w = pop(&BW_Successors[v]);
                
                if(v == s){          
                    if(numStartAtNode_temp[v] > 0){                        
                        aux = 
                           ((spEBW_ExcludeInit_temp[v*sizeSupra + w] + spEBW_Init_temp[v*sizeSupra + w])/numStartAtNode_temp[v])
                           *BW_delta[v]*inOutFluxRatios[v];
                    }else{
                        aux = 0;
                    }
                    spEBW_Init[v*sizeSupra + w] = spEBW_Init[v*sizeSupra + w] + aux;                    
                    BW_delta[w] = BW_delta[w] + aux;                                                            
                                        
                }else{                    
                    if(spNBW_Temp[v]>0){
                        aux = ((spEBW_ExcludeInit_temp[v*sizeSupra + w] + spEBW_Init_temp[v*sizeSupra + w])/spNBW_Temp[v])
                                    *(BW_delta[v]*BW_p[v])*inOutFluxRatios[v];
                    }else{
                        aux = 0;
                    }
                    spEBW_ExcludeInit[v*sizeSupra + w] = spEBW_ExcludeInit[v*sizeSupra + w] + aux;                    
                    BW_delta[w] = BW_delta[w] + aux;                    
                }                                
            }
                        
            if(v != s){
                spNBW[v] = spNBW[v] + BW_delta[v]*BW_p[v];
                numEndAtNode[v] = numEndAtNode[v] + BW_delta[v]*(1-BW_p[v]);                                       
            }            
        }
        // END: forward correction of the betweenness                        
    }
    
    //mxFree(spNBW_Temp);   
    //mxFree(spEBW_ExcludeInit_temp);   
    //mxFree(spEBW_Init_temp);     
    //mxFree(numEndAtNode_temp);      
    //mxFree(numStartAtNode_temp);                            
}

void fillDesinationsProbMatrix(
        double * possDestProbMatrix,
        const mxArray * possDest,
        const mxArray * possDestProb, int sizeSupra){
    
    int sourceIndex = -1;
    int targetIndex = -1;
    int numTargets = -1;
    mxArray * currPossDest = 0;
    double * currPossDestPr = 0;        
    double * currProbOfDestPr = 0;
    
    memset ( possDestProbMatrix, 0, sizeof(double)*sizeSupra*sizeSupra );
    
    for(sourceIndex = 0; sourceIndex<sizeSupra; sourceIndex++){
        currPossDest = mxGetCell(possDest, sourceIndex);
        
        numTargets = mxGetN(currPossDest);
        if(numTargets > 0 && mxGetM(currPossDest) > 0){ // check if the array is emtpy
            //get pointer to possible destinies
            currPossDestPr = mxGetPr(currPossDest);
            //get pointer to the probability matrix            
            currProbOfDestPr = mxGetPr(mxGetCell(possDestProb, sourceIndex));                        
            
            for(targetIndex = 0; targetIndex < numTargets ; targetIndex ++){    
                possDestProbMatrix[ ((int)currPossDestPr[targetIndex]-1)*sizeSupra + sourceIndex] = currProbOfDestPr[targetIndex];
            }
        }        
    }           
}

void recomputeJumpingPorbabilities(
        double * B_i_obs,
        double * e_i_obs,
        double * EB_Ext_obs,
        double * EB_Int_obs,
        double * pCondExtPacket, 
        double * probTransitExternPack,
        double * probTransitLocalPack,
        int sizeSupra){    
    
    int i;
    int j;
    
    // this goes for p
    for(i = 0 ; i < sizeSupra ; i++){
        if(B_i_obs[i]+e_i_obs[i] > 0){
            pCondExtPacket[i] = (B_i_obs[i])/(B_i_obs[i]+e_i_obs[i]);
        }else{
            pCondExtPacket[i] = 0;
        }
    }                 

    projectColumnsWithAddition(&M_totalNumRoutings[0], EB_Ext_obs, sizeSupra, sizeSupra);       
    for(j = 0 ; j < sizeSupra ; j++){            
        for(i = 0 ; i < sizeSupra ; i++){
            if(M_totalNumRoutings[i]>0){
                probTransitExternPack[i*sizeSupra + j] = EB_Ext_obs[i*sizeSupra + j]/M_totalNumRoutings[i];
            }
        }
    }
    
    projectColumnsWithAddition(&M_totalNumRoutings[0], EB_Int_obs, sizeSupra, sizeSupra);       
    for(j = 0 ; j < sizeSupra ; j++){            
        for(i = 0 ; i < sizeSupra ; i++){
            if(M_totalNumRoutings[i]>0){
                probTransitLocalPack[i*sizeSupra + j] = EB_Int_obs[i*sizeSupra + j]/M_totalNumRoutings[i];
            }
        }
    }    
    projectColumnsWithAddition(&cumProbTransitExternPack[0], probTransitExternPack, sizeSupra, sizeSupra);                
}

void projectRowsWithAddition(double * resultingVect, double * matrix, int numRows, int numColumns){
    int i,j;
    //copy the transient values in the Q numTransient
    memset(resultingVect,0,sizeof(double)*numColumns);
    for(i = 0; i < numRows ; i ++){        
        for(j = 0; j < numColumns ; j ++){
            resultingVect[j] = resultingVect[j] + matrix[i*numColumns + j];
        }
    }
}

void projectColumnsWithAddition(double * resultingVect, double * matrix, int numRows, int numColumns){
    int i,j;
    //copy the transient values in the Q numTransient
    memset(resultingVect,0,sizeof(double)*numColumns);        
    for(i = 0; i < numRows ; i ++){                
        for(j = 0; j < numColumns ; j ++){        
            resultingVect[i] = resultingVect[i] + matrix[i*numColumns + j];
        }
    }
}

void printMatrix(double * P, int rows, int cols){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < rows; rowIndex++){
        for(columnIndex = 0 ; columnIndex < cols; columnIndex++){            
            mexPrintf("%1.4f ",P[rowIndex*cols + columnIndex]);   
        }    
        mexPrintf("\n");   
    }        
}

void printMatrixWithIndex(double * P, int rows, int cols){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < rows; rowIndex++){
        for(columnIndex = 0 ; columnIndex < cols; columnIndex++){            
            mexPrintf("(%i,%i)%1.4f ",rowIndex,columnIndex,P[rowIndex*cols + columnIndex]);   
        }    
        mexPrintf("\n");   
    }        
}


void printMatrixInt(int * P, int rows, int cols){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < rows; rowIndex++){
        for(columnIndex = 0 ; columnIndex < cols; columnIndex++){            
            mexPrintf("%i ",P[rowIndex*cols + columnIndex]);   
        }    
        mexPrintf("\n");   
    }        
}

void printSquareMatrix(double * P, int sizeSupra){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < sizeSupra; rowIndex++){
        for(columnIndex = 0 ; columnIndex < sizeSupra; columnIndex++){            
            mexPrintf("%1.4f ",P[rowIndex*sizeSupra + columnIndex]);   
        }    
        mexPrintf("\n");   
    }        
}
