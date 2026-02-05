
// to compile use: mex -largeArrayDims SPDataPathDegeneration_staticMemory.c

#include "mex.h"
#include <math.h>
#include <string.h>
//#include "linkedList.h"
#define MAX_COLA 2000
#include "predQueueCircStatic.h"

void floydWarshall(const mxArray *G, double *D, staticQueue * P, double *PNumb, int numNodes, int numLayers);
void probOfExit(staticQueue * P, double * PNumb, double * costMat, double * numRoutingsLocalPacket, double * numPathsStart, double * numPathsEnd, int numNodes, int numLayers);
void localGenPath(staticQueue * P, double * PNumb, 
        double * numRoutingsLocalPacketCurr, double numPaths_ij, int I, int J, int source, int target, int numNodes, int numLayers);

/*
void printMatrix(double * P, int sizeSupra);
void printMatrix3D(double * P, int sizeSupra);
void localBetweenness(double *predMat, double * numPred, double * numPathMat, double * BWVect, double numPaths_ij, int I, int J, int currNode, int N);
void spBetweenness(double *predMat, double * numPred, double * numPathMat, double * BWVector, int sizeSupra);
void Path(double *predMat, double * numPred, double * numPathMat, double * BWVect, double numPaths_ij, int I, int J, int source, int target, int N);
*/
 
char outputString[1000];


/*
[Dist,NumPaths,Pred,PredProb] = SPDataPathDegeneration(sparse(A),5,1)
 * Dist: Pairwise distnaces -> column index source, row destination
 * NumPaths: Number of paths from a source to a destiny -> column index source, row destination
 * Pred: Predecessors matrix -> row index source, column destiny
 *          for the path that goes from source(row) to destiny (column) the previous node to reach the destiny 
 * PredProb: Cummulative probabilities of jumping to the a predecessor to have equal prob to go through any path
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{   
    
//     //output matrices
//     double *costMat,*PNumb;
//     queue ** P;
//     int numNodes = 0, numLayers = 0;
//     int sizeSupra = 0;
//     mwSignedIndex dimsP[3];
//     
//     // get and check inputs
//     if (nrhs != 3) mexErrMsgTxt( "Only 3 input arguments allowed." );
//     if (nlhs != 3) mexErrMsgTxt( "Only 3 output argument allowed." );
//     sizeSupra = mxGetN( prhs[0] );
//     if (mxGetM( prhs[0] ) != sizeSupra) mexErrMsgTxt( "Input matrix G needs to be square." );
//     if(mxIsSparse(prhs[0])==0) mexErrMsgTxt( "Distance Matrix must be sparse" );
// 
//     numNodes = mxGetPr(prhs[1])[0];
//     numLayers = mxGetPr(prhs[2])[0];   
//     
//     if((numNodes*numLayers)!=sizeSupra) mexErrMsgTxt( "Incorrect size of the input matrix" );
//     
//     costMat = mxGetPr(mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL ));
//     PNumb = mxGetPr(mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL ));
//     P = (queue **)mxMalloc(sizeof(queue *)*sizeSupra*sizeSupra);       
//     
//     mexEvalString("fprintf('Start betweenness algorithm\\n');");         
//     floydWarshall(prhs[0], costMat, P, PNumb, numNodes, numLayers);
//  
//     //mexEvalString("fprintf('Number of paths\\n');");             
//     //printMatrix(PNumb, sizeSupra);    
//     //mexEvalString("fprintf('Cost matrix\\n');");             
//     //printMatrix(costMat, sizeSupra);    
//     //printMatrix3D(P, sizeSupra);
//     //printMatrix(PElems, sizeSupra);    
//     
//     int i;
//     double * numRoutingsLocalPacket; 
//     double * numPathsStart; 
//     double * numPathsEnd; 
//     //double * numOutPaths;
//     //double * probOfCont;
//     //double * totalPacketsPassNode;
//     
//     plhs[0] = mxCreateDoubleMatrix( numNodes*numLayers, numNodes*numLayers, mxREAL );
//     plhs[1] = mxCreateDoubleMatrix( numNodes*numLayers, 1, mxREAL );
//     plhs[2] = mxCreateDoubleMatrix( numNodes*numLayers, 1, mxREAL );
// 
//     numRoutingsLocalPacket = mxGetPr(plhs[0]);
//     numPathsStart = mxGetPr(plhs[1]);
//     numPathsEnd = mxGetPr(plhs[2]);
//     
//     probOfExit(P,PNumb,costMat,numRoutingsLocalPacket,numPathsStart,numPathsEnd,numNodes,numLayers);
//     mexEvalString("fprintf('The end!\\n');");       
    
    //output matrices
    mxArray * curArray;
    double * curArrayPtr;
    mxArray * curArrayProb;
    double * curArrayProbPtr;    
    double *costMat,*PNumb;
    double numPaths;
    staticQueue * P;
    int numNodes = 0, numLayers = 0;
    int sizeSupra = 0, i, j, c, numEls, z;
    mwSignedIndex dimsP[3];
    mxArray * POptionsCell;
    mxArray * POptionsProb;
    
    // get and check inputs
    if (nrhs != 3) mexErrMsgTxt( "Only 3 input arguments allowed." );
    if (nlhs > 10) mexErrMsgTxt( "Only 10 output argument allowed." );
    sizeSupra = mxGetN( prhs[0] );
    if (mxGetM( prhs[0] ) != sizeSupra) mexErrMsgTxt( "Input matrix G needs to be square." );
    if(mxIsSparse(prhs[0])==0) mexErrMsgTxt( "Distance Matrix must be sparse" );

    numNodes = mxGetPr(prhs[1])[0];
    numLayers = mxGetPr(prhs[2])[0];   
    
    if((numNodes*numLayers)!=sizeSupra) mexErrMsgTxt( "Incorrect size of the input matrix" );        
    
    plhs[0] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );
    costMat = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );
    PNumb = mxGetPr(plhs[1]);
    P = (staticQueue *)mxMalloc(sizeof(staticQueue)*sizeSupra*sizeSupra);
    
    
    mexEvalString("fprintf('Start betweenness algorithm\\n');");       
    floydWarshall(prhs[0], costMat, P, PNumb, numNodes, numLayers);    
    mexEvalString("fprintf('End betweenness algorithm\\n');");   
    
    
    /****************************************************************************/
    /********************** Output data for the theoretical *********************/
    /****************************************************************************/
//     double * numRoutingsLocalPacket; 
//     double * numPathsStart; 
//     double * numPathsEnd; 
//     
//     plhs[7] = mxCreateDoubleMatrix( numNodes*numLayers, numNodes*numLayers, mxREAL );
//     plhs[8] = mxCreateDoubleMatrix( numNodes*numLayers, 1, mxREAL );
//     plhs[9] = mxCreateDoubleMatrix( numNodes*numLayers, 1, mxREAL );
// 
//     numRoutingsLocalPacket = mxGetPr(plhs[7]);
//     numPathsStart = mxGetPr(plhs[8]);
//     numPathsEnd = mxGetPr(plhs[9]);
//     
//     probOfExit(P,PNumb,costMat,numRoutingsLocalPacket,numPathsStart,numPathsEnd,numNodes,numLayers);
//     mexEvalString("fprintf('The end!\\n');");             
        
    /****************************************************************************/
    /********************** Output data for the simulations *********************/
    /****************************************************************************/
    
    //copy predecessor list into a cell array
    POptionsCell = mxCreateCellMatrix(sizeSupra,sizeSupra);
    POptionsProb = mxCreateCellMatrix(sizeSupra,sizeSupra);
    plhs[2] = POptionsCell;
    plhs[3] = POptionsProb;
    for(i=0;i<sizeSupra;i++){
        for(j=0;j<sizeSupra;j++){            
            //sprintf(outputString,"fprintf('i=%i,j=%i\\n')",i,j);
            //mexEvalString(outputString);                  
            
            curArray = mxCreateDoubleMatrix( numElems(&P[i*sizeSupra + j]),1, mxREAL );
            curArrayProb = mxCreateDoubleMatrix( numElems(&P[i*sizeSupra + j]),1, mxREAL );
            curArrayPtr = mxGetPr(curArray);
            curArrayProbPtr = mxGetPr(curArrayProb);
            
            numEls = numElems(&P[i*sizeSupra + j]);            
            numPaths = 0; //in order to create a probability of jumping
            for(c = 0; c < numEls ; c++){
                z = get(&P[i*sizeSupra + j]);
                curArrayPtr[c] = (double)z+1;
                
                numPaths = numPaths + PNumb[z*sizeSupra + i];
                curArrayProbPtr[c] = PNumb[z*sizeSupra + i];
            }
    
            if(numEls > 0){
                if(numPaths > 0){
                    curArrayProbPtr[0] = curArrayProbPtr[0]/numPaths;
                    for(c = 1; c < numEls ; c++){
                        curArrayProbPtr[c] = curArrayProbPtr[c-1] + curArrayProbPtr[c]/numPaths;
                        //curArrayProbPtr[c] = curArrayProbPtr[c]/numPaths;
                    }
                }else{
                    if(curArrayPtr[0]==(i+1)){
                        curArrayProbPtr[0] = 1;
                    }                    
                }
            }
        
            mxSetCell(POptionsProb,j*sizeSupra + i,curArrayProb);            
            mxSetCell(POptionsCell,j*sizeSupra + i,curArray);
        }    
    }
    mexEvalString("fprintf('End Predecessor list\\n');");   
    
    //generate structures to choose the best source and destiny considering the different layers
    int s,d,l1,l2;
    mxArray * curSrc, * curDst;
    double * curSrcPr, * curDstPr;
    mxArray * curPathProbCum;
    double * curPathProbCumPr;
    double curCost;
    double tolerance = 0.0000000001; 
    staticQueue sources;
    staticQueue destinies;            
    staticQueue numPathsS_DList;
    int totalPaths = 0;
    double minCost = mxGetInf();
    mxArray * spSrcsMult = mxCreateCellMatrix(sizeSupra,sizeSupra);
    plhs[4] = spSrcsMult; 
    mxArray * spDestMult = mxCreateCellMatrix(sizeSupra,sizeSupra);
    plhs[5] = spDestMult; 
    mxArray * spPathProbCum = mxCreateCellMatrix(sizeSupra,sizeSupra);
    plhs[6] = spPathProbCum;
    
    int z1, z2, p;
    int sM,dM;
    for(s = 0; s<sizeSupra; s++){
        for(d = 0; d<sizeSupra; d++){
            //if(s != d){
                //create matrix to choose the best source and destination nodes
                //initialize the queues
                create(&sources);
                create(&destinies);            
                create(&numPathsS_DList);
                minCost = mxGetInf();
                totalPaths = 0;
                
                for(l1 = 0; l1<numLayers; l1++){
                    for(l2 = 0; l2<numLayers; l2++){

                        sM = (l1*numNodes + (s%numNodes));
                        dM = (l2*numNodes + (d%numNodes));                    

                        curCost = costMat[sM*sizeSupra + dM];                    

                        //sprintf(outputString,"fprintf('doing source = %i, new destiny = %i, cost = %f\\n')",sM,dM,curCost);
                        //mexEvalString(outputString);                         


                        if(minCost>curCost && (fabs(minCost - curCost) > tolerance)){
                            create(&sources);
                            create(&destinies);
                            create(&numPathsS_DList);
                            //add the new minimum
                            add(&sources, sM); 
                            add(&destinies,dM);                                                                       
                            add(&numPathsS_DList,(int)PNumb[sM*sizeSupra + dM]); 
                            totalPaths = PNumb[sM*sizeSupra + dM];
                            minCost = curCost;       

                            //sprintf(outputString,"fprintf('\\tnew source = %i, new destiny = %i, cost = %f, totalPaths = %i\\n')",sM,dM,curCost,totalPaths);
                            //mexEvalString(outputString);                         


                        }else if(fabs(minCost - curCost) < tolerance){
                            add(&sources, sM);
                            add(&destinies,dM);                        
                            add(&numPathsS_DList,PNumb[sM*sizeSupra + dM]); 
                            totalPaths = totalPaths + PNumb[sM*sizeSupra + dM];                        

                            //sprintf(outputString,"fprintf('\\tnew path source = %i, new destiny = %i, cost = %f\\n')",sM,dM,curCost);
                            //mexEvalString(outputString);                         
                        }                    
                    }
                }

                //create sources array
                numEls = numElems(&sources);            
                curSrc = mxCreateDoubleMatrix( numEls,1, mxREAL );    
                curSrcPr = mxGetPr(curSrc);
                curDst = mxCreateDoubleMatrix( numEls,1, mxREAL );   
                curDstPr = mxGetPr(curDst);
                curPathProbCum = mxCreateDoubleMatrix( numEls,1, mxREAL );
                curPathProbCumPr = mxGetPr(curPathProbCum);
                for(c = 0; c < numEls ; c++){
                    z1 = get(&sources);
                    z2 = get(&destinies);
                    p = get(&numPathsS_DList);                                       
                    curSrcPr[c] = (double)z1+1;
                    curDstPr[c] = (double)z2+1;
                    
                    if(c == 0){
                        curPathProbCumPr[0] = (double)p / (double)totalPaths;
                        //sprintf(outputString,"fprintf('\\tcurPathProbCumPr[0] = %f, p = %i, totalPaths = %i\\n')",curPathProbCumPr[0],p,totalPaths);
                        //mexEvalString(outputString);                                                 
                    }else{
                        curPathProbCumPr[c] = curPathProbCumPr[c-1] + (double)p / (double)totalPaths;
                        //sprintf(outputString,"fprintf('\\tcurPathProbCumPr[c] = %f\\n')",curPathProbCumPr[c]);
                        //mexEvalString(outputString);                         
                    }
                }            
                mxSetCell(spSrcsMult,s*sizeSupra + d,curSrc); 
                mxSetCell(spDestMult,s*sizeSupra + d,curDst); 
                mxSetCell(spPathProbCum,s*sizeSupra + d,curPathProbCum); 
                
            //}
        }
    }
    mexEvalString("fprintf('End all list\\n');");       
    
}

void floydWarshall(const mxArray *G, double *D, staticQueue * P, double *PNumb, int numNodes, int numLayers){
    double INF=mxGetInf(),SMALL=mxGetEps();
    double tolerance = 0.0000000001;       
    int sizeSupra = mxGetN(G);
    int rowIndex = 0, columnIndex = 0, i = 0, numEls = 0, c = 0;
    int newPred;
    double *Gpr = mxGetPr(G);
    size_t *Gir = mxGetIr(G);
    size_t *Gjc = mxGetJc(G);    
    
    mexEvalString("fprintf('Start of floyd-warshall!\\n');");         
    
    //initialize everything to initial value
    for(rowIndex = 0 ; rowIndex < sizeSupra; rowIndex++){
        
        //sprintf(outputString,"fprintf('\\tInitialization iteration %i of %i!\\n');",rowIndex,sizeSupra);
        //mexEvalString(outputString);                 
        
        for(columnIndex = 0 ; columnIndex < sizeSupra; columnIndex++){
            //initial dists to inf
            D[columnIndex*sizeSupra + rowIndex] = INF;
            //number of paths bw pairs = 0
            PNumb[columnIndex*sizeSupra + rowIndex] = 0;            
            create(&P[rowIndex*sizeSupra + columnIndex]);            

            //This would be equivalent to having a self loop
            if(columnIndex == rowIndex){
                D[columnIndex*sizeSupra + rowIndex] = 0;
                PNumb[columnIndex*sizeSupra + rowIndex] = 1;            
                add(&P[rowIndex*sizeSupra + columnIndex], columnIndex);
            }            
        }    
    }   
    
    //initialize the entries that contain a value
    for(columnIndex = 0; columnIndex<sizeSupra; columnIndex++){                         
        //mexPrintf("closest = %i, startInd = %i, endInd = %i\n",closest, startInd,endInd);
        //sprintf(outputString,"fprintf('\\tLoading network data %i of %i!\\n');",columnIndex,sizeSupra);
        //mexEvalString(outputString);                 
        
        if( Gjc[columnIndex]!=Gjc[columnIndex+1] )
            for( i=Gjc[columnIndex]; i<=(Gjc[columnIndex+1]-1); i++ ) {
                rowIndex = Gir[i];
                D[columnIndex*sizeSupra + rowIndex] = Gpr[i];
                PNumb[columnIndex*sizeSupra + rowIndex] = 1; 
                add(&P[rowIndex*sizeSupra + columnIndex], rowIndex);
            }    
    }        
    
    mexEvalString("fprintf('Initialization of floyd-warshall finished!\\n');");         
    
    int K = 0;
    int I = 0;
    int J = 0;
    int N = sizeSupra;
    //the algorithm
    for(K=0; K<N; K++){        
        sprintf(outputString,"fprintf('\\tK=%i!\\n');",K);
        mexEvalString(outputString);         

        for(I=0; I<N ; I++){
            if(I != K){
                for(J=0 ; J<N ; J++){
                    if(J!=K && J!=I){   
                        if(D[K*N+I] < INF && D[J*N+K] < INF){
                            if(D[K*N+I] + D[J*N+K] < D[J*N+I]){
                                D[J*N+I] = D[K*N+I] + D[J*N+K];
                                PNumb[J*N+I] = PNumb[K*N+I]*PNumb[J*N+K];
                                
                                //deal with the list of paths --> Preds(I, J) := Preds(K, J);
                                create(&P[I*N + J]);                                
                                numEls = numElems(&P[K*sizeSupra + J]);
                               
                                for(c=0;c<numEls;c++){
                                    add(&P[I*N + J], consultIndex(&P[K*sizeSupra + J], c)); //ja que primer s'ha de passar per k
                                }
                                
                            }else if(fabs(D[K*N+I] + D[J*N+K] - D[J*N+I])<tolerance){
                                PNumb[J*N+I] = PNumb[J*N+I] + PNumb[K*N+I]*PNumb[J*N+K];
                                
                                /*
                                sprintf(outputString,"fprintf('New path from %i to %i that goes through = %i found\\n')"
                                        ,I,J,K);
                                mexEvalString(outputString);
                                */
                                //copy the predecessors of K
                                numEls = numElems(&P[K*sizeSupra + J]);                               
                                for(c=0;c<numEls;c++){                   
                                    newPred = consultIndex(&P[K*sizeSupra + J], c);
                                    if(exist(&P[I*N + J], newPred) == FALSE){
                                        add(&P[I*N + J], newPred);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void probOfExit(staticQueue * P, double * PNumb, double * costMat, double * numRoutingsLocalPacket, double * numPathsStart, double * numPathsEnd, int numNodes, int numLayers){        
    int source = 2, currSourceNode = 0, currTargetNode = 0;
    int target = 0;
    int i,j;
    int sizeSupra = numNodes*numLayers;
    int minCostSourceNode[numLayers*numLayers];    
    int minCostSourceLayer[numLayers*numLayers];        
    int minCostTargetNode[numLayers*numLayers];    
    int minCostTargetLayer[numLayers*numLayers];  
    int numMins = 0;
    double minCost = 0;
    double currCost = -1;
    int pathIndex = 0, numberOfPathsInTotal = 0;
    int layerSce = 0, layerTgt = 0, row = 0, column = 0;
    double INF = mxGetInf();
    double * numRoutingsLocalPacketCurr = mxGetPr(mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL ));
        
    //double * localSPRouting = mxGetPr(mxCreateDoubleMatrix( numNodes*numLayers, numNodes*numLayers, mxREAL ));
        
    memset( (void *)numRoutingsLocalPacket, 0, sizeof(double)*numNodes*numLayers*numNodes*numLayers );
    memset( (void *)numPathsStart, 0, sizeof(double)*numNodes*numLayers );
    memset( (void *)numPathsEnd, 0, sizeof(double)*numNodes*numLayers );
    
    for(source = 0; source < numNodes ; source ++){
        //source = 0;
        for(target = 0; target < numNodes ; target ++){   
            //source = 4;
            //target = 203;
            if(source != target){

                //sprintf(outputString,"fprintf('->source = %i, target = %i\\n');",source,target);
                //mexEvalString(outputString);                                 
                
                minCost = INF;
                numMins = 0;
                
                for(layerSce = 0; layerSce < numLayers ; layerSce ++){
                    for(layerTgt = 0; layerTgt < numLayers ; layerTgt ++){
                        
                        row = (layerSce*numNodes+source);
                        column = (layerTgt*numNodes+target);     
                        currCost = costMat[row*sizeSupra + column];

                        if(currCost < INF){
                            if(currCost < minCost){
                                minCostSourceNode[0] = source;    
                                minCostSourceLayer[0] = layerSce;        
                                minCostTargetNode[0] = target;    
                                minCostTargetLayer[0] = layerTgt;                             
                                numMins = 1;        
                                minCost = costMat[row*sizeSupra + column];
                            }else if(currCost == minCost){
                                minCostSourceNode[numMins] = source;    
                                minCostSourceLayer[numMins] = layerSce;        
                                minCostTargetNode[numMins] = target;    
                                minCostTargetLayer[numMins] = layerTgt;                             
                                numMins ++;
                            }                          
                        }
                    }
                }                                        
        
                numberOfPathsInTotal = 0;
                for(pathIndex = 0; pathIndex < numMins; pathIndex++){
                    currSourceNode = minCostSourceLayer[pathIndex]*numNodes + minCostSourceNode[pathIndex];
                    currTargetNode = minCostTargetLayer[pathIndex]*numNodes + minCostTargetNode[pathIndex];
                    numberOfPathsInTotal += (int)PNumb[currSourceNode*sizeSupra + currTargetNode];
                }
                                            
                
                //sprintf(outputString,"fprintf('source = %i, target = %i, numPahts = %i\\n');",source,target,numMins);
                //mexEvalString(outputString);                                 
                
                for(pathIndex = 0; pathIndex < numMins; pathIndex++){

                    currSourceNode = minCostSourceLayer[pathIndex]*numNodes + minCostSourceNode[pathIndex];
                    currTargetNode = minCostTargetLayer[pathIndex]*numNodes + minCostTargetNode[pathIndex];
                    
                    //mexPrintf("path from %i to %i, cost = %f ,normalization %i/%i\n", currSourceNode,currTargetNode,minCost,
                    //        (int)PNumb[currSourceNode*sizeSupra + currTargetNode],numberOfPathsInTotal);                    
                    
                    //sprintf(outputString,"fprintf('\\t%i\\n');",(int)currSourceNode);
                    //mexEvalString(outputString);   
                    
                    memset ( (void *)numRoutingsLocalPacketCurr, 0, sizeof(double)*numNodes*numLayers*numNodes*numLayers );                        
                    
//                      if(source == 48){
//                          mexPrintf("path from %i to %i, cost = %f ,normalization %i/%i\n", currSourceNode + 1,currTargetNode + 1,minCost,
//                                  (int)PNumb[currSourceNode*sizeSupra + currTargetNode],numberOfPathsInTotal);                    
//                      }
                        
                    numPathsStart[currSourceNode] = numPathsStart[currSourceNode] + 
                            (PNumb[currSourceNode*sizeSupra + currTargetNode]/numberOfPathsInTotal);
                    numPathsEnd[currTargetNode] = numPathsEnd[currTargetNode] + 
                            (PNumb[currSourceNode*sizeSupra + currTargetNode]/numberOfPathsInTotal);                    

                    
                    
                    localGenPath(P,PNumb,numRoutingsLocalPacketCurr,
                        PNumb[currSourceNode*sizeSupra + currTargetNode],
                        currSourceNode,
                        currTargetNode,
                        currSourceNode,
                        currTargetNode,
                        numNodes,numLayers);    
                    
                    for(i=0;i<sizeSupra;i++){
                        for(j=0;j<sizeSupra;j++){
                            numRoutingsLocalPacket[i*sizeSupra + j] = numRoutingsLocalPacket[i*sizeSupra + j] +
                                    (PNumb[currSourceNode*sizeSupra + currTargetNode]/numberOfPathsInTotal)*numRoutingsLocalPacketCurr[i*sizeSupra + j];
                            
                        }                    
                    }                                       
                }
            }
            //return;
        }
        //return;
    }
}

void localGenPath(staticQueue * P, double * PNumb, 
        double * numRoutingsLocalPacketCurr, double numPaths_ij, int I, int J, int source, int target, int numNodes, int numLayers){
    
    int intermediate;
    int intermediateIndex;
    int sizeSupra = numNodes*numLayers;
    int numElements;
    
    //number of elements in the list
    numElements = numElems(&P[I*sizeSupra + J]);    
    //mexPrintf("\t numIntermediate(%i,%i) %i\n",I,J,(int)numElements);    
    
    if(numElements > 0){
        intermediate = consultIndex(&P[I*sizeSupra + J], 0);
        
        if(numElements == 1 && intermediate == source){               
            //sprintf(outputString,"fprintf('\\t %i - %i\\n');",I,J);
            //mexEvalString(outputString);             
            if(J != target){
                numRoutingsLocalPacketCurr[sizeSupra*I + J] = PNumb[J*sizeSupra + target]/PNumb[I*sizeSupra + target];
            }

            if(J == target){
                numRoutingsLocalPacketCurr[sizeSupra*I + J] = 1;
            }
        }
                
            
        if(intermediate != J && intermediate != I){
            for(intermediateIndex = 0; intermediateIndex < numElements; intermediateIndex ++){

                intermediate = consultIndex(&P[I*sizeSupra + J], intermediateIndex);

                //extend the search
                localGenPath(P,PNumb,numRoutingsLocalPacketCurr,numPaths_ij,I,intermediate,source,target,numNodes,numLayers);                                
                //mexPrintf("\t\t %i\n",(int)intermediate);        
                localGenPath(P,PNumb,numRoutingsLocalPacketCurr,numPaths_ij,intermediate,J,source,target,numNodes,numLayers);                        
            }
        }
        
    }
}



void printMatrix(double * P, int sizeSupra){
    int rowIndex = 0, columnIndex = 0;
    
    //print the matrix
    for(rowIndex = 0 ; rowIndex < sizeSupra; rowIndex++){
        for(columnIndex = 0 ; columnIndex < sizeSupra; columnIndex++){            
            mexPrintf("%f ",P[rowIndex*sizeSupra + columnIndex]);   
        }    
        mexPrintf("\n");   
    }    
    
}


void printMatrix3D(double * P, int sizeSupra){
    int rowIndex = 0, columnIndex = 0, deepIndex = 0;
    
    //print the matrix
    for(deepIndex = 0; deepIndex < sizeSupra; deepIndex ++){
        mexPrintf("deepIndex = %i\n\t", deepIndex);
        for(rowIndex = 0 ; rowIndex < sizeSupra; rowIndex++){
            for(columnIndex = 0 ; columnIndex < sizeSupra; columnIndex++){            
                mexPrintf("%f ",P[rowIndex*sizeSupra*sizeSupra + columnIndex*sizeSupra + deepIndex]);   
            }    
            if(rowIndex+1 < sizeSupra)
                mexPrintf("\n\t");   
            else
                mexPrintf("\n");   
        }    
    }
    
}

// Distances.all := (others => (others => Infinity));
//   Paths.all := (others => (others => 0));
//   for I in 1..N loop
//     if Pgr(I).From /= null then
//       for K in Pgr(I).From'Range loop
//         J := Pgr(I).From(K).Index;
//         if I /= J then
//           Distances(I, J) := 1;
//           Paths(I, J) := 1;
//         end if;
//       end loop;
//     end if;
//     Distances(I, I) := 0;
//   end loop;
// 
//   for K in 1..N loop
//     for I in 1..N loop
//       if I /= K then
//         for J in 1..N loop
//           if J /= K and J /= I then
//             if Distances(I, K) < Infinity and Distances(K, J) < Infinity then
//               if Distances(I, K) + Distances(K, J) < Distances(I, J) then
//                 Distances(I, J) := Distances(I, K) + Distances(K, J);
//                 Paths(I, J) := Paths(I, K) * Paths(K, J);
//               elsif Distances(I, K) + Distances(K, J) = Distances(I, J) then
//                 Paths(I, J) := Paths(I, J) + Paths(I, K) * Paths(K, J);
//               end if;
//             end if;
//           end if;
//         end loop;
//       end if;
//     end loop;
//   end loop;
//   


//  READ SPARSE MATRIX EXAMPLE
//     double *Gpr = mxGetPr(prhs[0]);
//     size_t *Gir = mxGetIr(prhs[0]);
//     size_t *Gjc = mxGetJc(prhs[0]);    
//     
//     int columnIndex = 0;
//     int i = 0;
//     long int whichNeigh = 0;
//     double arcLength = 0;    
//     
//     for(columnIndex = 0; columnIndex<sizeSupra; columnIndex++){
//         long int startInd = Gjc[ columnIndex   ];
//         long int endInd   = Gjc[ columnIndex+1 ] - 1;
//                          
//         //mexPrintf("closest = %i, startInd = %i, endInd = %i\n",closest, startInd,endInd);
//         for( i=startInd; i<=endInd; i++ ) {
//             whichNeigh = Gir[ i ];  //neigbour of closest ?
//             arcLength = Gpr[ i ];   
//             mexPrintf("(%i,%i) -> %f\n",whichNeigh+1,columnIndex+1,arcLength);         
//         }    
//     }