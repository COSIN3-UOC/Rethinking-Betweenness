//to compile: mex -largeArrayDims cSPCongestion_statMem_dir_weighted_local.c

#include "mex.h"
#define MAX_COLA 501000
#include "packetQueueCircStatic.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifndef __RANDOM__
#define __RANDOM__
    // Random generators Libs
    #include "./RandomLib/rvgs.c"
    #include "./RandomLib/rngs.c"
#endif

        
#define TRUE 1
#define FALSE 0
#define DEBUG 0
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))       
               
char outputString[1000];
    
double FastGetNextRandomNode_Linear(double * cumArray, long length, double r){
    int k = 0;

    if(length > 0){
        if(cumArray[length-1] > 0){
            for(k=0; k<length; k++){
                //sprintf(outputString,"fprintf('%i=%f,')",k,SupraTransitionCumMatrix[k*(numNodes*numLayers) + supraNode]);
                //mexEvalString(outputString);                  
                if(cumArray[k]>=r){
                    return k;
                }
            }        
            //sprintf(outputString,"fprintf('\\n')",k,SupraTransitionCumMatrix[k*(numNodes*numLayers) + supraNode]);
            //mexEvalString(outputString);                  
        }else{
            sprintf(outputString,"fprintf('FastGetNextRandomNode_Linear: all zeros in array\\n')");
            mexEvalString(outputString);          
            return -1;
        }
    }
    return -1;
}

//c = cSPCongestion_statMem(0.2,1,Dist,NumPaths,Pred,4,1,timeSteps,seed);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int timeSteps, seed;
    int * genRate;
    double * genRateProb;
    int * delivRate;
    double * delivRatePtr;
    staticQueue * networkQueues;        
    int p,n,t,l,supraNodeId;
    int numNodes, numLayers;
    packet newPacket;
    packet * newPacket_ptr;
    int delivered, conToDest;
    double r = -1;
    int nextDest;
    double * numPackets;
    const mxArray * srcDestProbMult;
    double probOfJump;
    int i,j;
    int randIndex = -1;
    int randDest = -1;
    int reqDataPack = 0;
    int sizeSupra = 0;
    mxArray * curCumProb;
    double * curCumProbPr;
    double * curCumPr;
    size_t numPred;
    double tolerance = 0.0000000001;
    double * currAccessibleDestPr;
    double * currProbOfDestPr;
    const mxArray * currAccessibleDest;
    const mxArray * currProbOfDest;
    long int deliveredPacketsTotal = 0;
    long int createdPackets = 0;
    int edgeSource = -1, edgeDest = -1;    
    int packetType=0;
    
    if ((nrhs != 15) || (nlhs > 12)) {
        mexPrintf("Error! Expecting exactly 16 input and up to 12 output argument!\n");
        return;
    }           
    
    //read input and generate output vars
    numNodes = (int)mxGetPr(prhs[11])[0];
    numLayers = (int)mxGetPr(prhs[12])[0];     
    timeSteps = (int)mxGetPr(prhs[13])[0];
    seed = mxGetPr(prhs[14])[0];        
    sizeSupra = numNodes*numLayers;
    // read the generation and deliver rates
    double * genRateOrig = mxGetPr(prhs[0]);
    genRateProb = (double *)mxMalloc(sizeof(double)*sizeSupra);
    genRate = (int *)mxMalloc(sizeof(int)*sizeSupra);
    delivRate = (int *)mxMalloc(sizeof(int)*sizeSupra);
    delivRatePtr = mxGetPr(prhs[1]);
    for(i = 0; i<numNodes; i++){
        genRate[i] = (int)genRateOrig[i];
        genRateProb[i] = genRateOrig[i] - (double)genRate[i];
        delivRate[i] = (int)delivRatePtr[i];
    }      
    if(DEBUG){
        for(i = 0; i<numNodes; i++){
            if((genRateProb[i] + (double)genRate[i]) > tolerance){
                mexPrintf("val %i = (%i,%2.2f)!\n",i,delivRate[i],genRateProb[i] + genRate[i]);
            }
        }         
    }
    
    //double * Dist = mxGetPr(prhs[2]);
    //double * numPaths = mxGetPr(prhs[3]);    
    const mxArray * Pred = prhs[4];    //això farà falta canviar-ho potser
    const mxArray * PredProb = prhs[5];   // això també     
    const mxArray * srcMult = prhs[6];    
    const mxArray * dstMult = prhs[7];    
    const mxArray * probMult = prhs[8];                    
    const mxArray * possDest = prhs[9];    
    const mxArray * possDestCumProb = prhs[10];    
    
    PlantSeeds(seed);
    //output variables
    plhs[0] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    numPackets = mxGetPr(plhs[0]);   
    memset ( numPackets, 0, sizeof(double)*timeSteps*sizeSupra); 

    plhs[1] = mxCreateDoubleMatrix(1, timeSteps, mxREAL);    
    double * numPacketOutUnitTime = mxGetPr(plhs[1]);   
    memset ( numPacketOutUnitTime, 0, sizeof(double)*timeSteps ); 

    plhs[2] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    double * B_i_ext = mxGetPr(plhs[2]);   
    memset ( B_i_ext, 0, sizeof(double)*sizeSupra*timeSteps);         

    plhs[3] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    double * e_i_ext = mxGetPr(plhs[3]);   
    memset ( e_i_ext, 0, sizeof(double)*sizeSupra*timeSteps);         

    plhs[4] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    double * B_i = mxGetPr(plhs[4]);   
    memset ( B_i, 0, sizeof(double)*sizeSupra*timeSteps);        

    plhs[5] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    double * e_i = mxGetPr(plhs[5]);   
    memset ( e_i, 0, sizeof(double)*sizeSupra*timeSteps);             

    plhs[6] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
    double * s_i = mxGetPr(plhs[6]);   
    memset ( s_i, 0, sizeof(double)*sizeSupra*timeSteps);             
       
    plhs[7] = mxCreateDoubleMatrix(sizeSupra, sizeSupra, mxREAL);    
    double * PExp = mxGetPr(plhs[7]);   
    memset ( PExp, 0, sizeof(double)*sizeSupra*sizeSupra);                 
 
     plhs[8] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
     double * sigma = mxGetPr(plhs[8]);   
     memset ( sigma, 0, sizeof(double)*sizeSupra*timeSteps);       
 
     plhs[9] = mxCreateDoubleMatrix(sizeSupra, timeSteps, mxREAL);    
     double * extractedPacketsTemp = mxGetPr(plhs[9]);   
     memset ( extractedPacketsTemp, 0, sizeof(double)*sizeSupra*timeSteps);            

     
//     //vars to debug
     plhs[10] = mxCreateDoubleMatrix(sizeSupra,1, mxREAL);    
     plhs[11] = mxCreateDoubleMatrix(sizeSupra,1, mxREAL);    
     double * extractedPackets = mxGetPr(plhs[10]);
     double * packetsThatContinue = mxGetPr(plhs[11]);
     memset ( extractedPackets, 0, sizeof(double)*sizeSupra ); 
     memset ( packetsThatContinue, 0, sizeof(double)*sizeSupra ); 
    
    
    //test parameter sizes
    //sizeSupra = mxGetN( prhs[2] );
    //if((numNodes*numLayers)!=sizeSupra) mexErrMsgTxt( "Incorrect size of the input matrix" );    
        
    t = 0;           
    
    //vars of the simulation
    networkQueues = mxMalloc(sizeSupra*sizeof(staticQueue));    
    for(n=0; n<sizeSupra ; n++){
        create(&networkQueues[n]);
    }                   
        
    numPackets[0] = 0;
    //count packets in the network
    for(n=0; n<sizeSupra; n++){
        numPackets[0] = numPackets[0] + numElems(&networkQueues[n]);
    }
    if(DEBUG){
        sprintf(outputString,"fprintf('Num packets @%i in system %i\\n')",0,(int)numPackets[0]);
        mexEvalString(outputString);           
    }   
    

    for(t=1;t<timeSteps;t++){
        
        if(t%10000 == 0){
            sprintf(outputString,"fprintf('Time %i of %i\\n')",t,timeSteps);
            mexEvalString(outputString);           
        }
        
        //if(t < 10){
            //generate packets
            //totalPackets[t] = 0;
            for(n=0; n<sizeSupra; n++){
                //if genRate is lower than one, we interpret that
                //is a probability, otherwise is number of packets            
                reqDataPack = 0;
                if(genRateProb > 0){
                    r = Uniform(0,1);
                    if(r<=genRateProb[n]){
                        reqDataPack = 1;
                    }else{
                        reqDataPack = 0;
                    }
                }

                reqDataPack = genRate[n] + reqDataPack;

                //generate the required data packages
                for(p=0; p<reqDataPack; p++){     
                    
                    //get the possible destinations
                    currAccessibleDest = mxGetCell(possDest, n);
                    
                    //it is possible to go somewhere                   
                    if(mxGetN(currAccessibleDest) > 0 && mxGetM(currAccessibleDest) > 0){
                                                
                        currAccessibleDestPr = mxGetPr(currAccessibleDest);                        
                        
                        //choose the destination according to the probabilities
                        currProbOfDest = mxGetCell(possDestCumProb, n);
                        currProbOfDestPr = mxGetPr(currProbOfDest);                        
                        
                        r = Uniform(0,1);     
                        randIndex = (int) FastGetNextRandomNode_Linear(
                                currProbOfDestPr,
                                mxGetN(currProbOfDest),
                                r);                         
                        
                        randDest = (int)currAccessibleDestPr[randIndex]-1;                                               
                        
                        if(DEBUG){
                            sprintf(outputString,"fprintf('r = %f; source %i has chosen dest %i\\n')",
                                    r,n,(int)randDest+1);
                            mexEvalString(outputString);                         
                        }
                            
                        if(randIndex != -1){
                            memset(&newPacket,0,sizeof(packet));

                            newPacket.source = (int)n;
                            newPacket.destination = (int)randDest;                   

                            //newPacket.source = 1;
                            //newPacket.destination = 1;                   

                            newPacket.maxSteps = -1;                         
                            //newPacket.maxSteps = Dist[newPacket.source*sizeSupra + newPacket.destination];                                    
                            newPacket.whereAmI = newPacket.source;                
                            newPacket.previous = newPacket.whereAmI;                           
                            //newPacket.deliveredOnTime = t-1;                                                                       
                            newPacket.deliveredOnTime = t;                                                                       

                           //effective betweenness
                           //if(nextDest != newPacket.destination){
                           //    B_i[t*sizeSupra + newPacket.source] ++; 
                           //}else{
                           //    e_i[t*sizeSupra + newPacket.source] ++; 
                           //}                   

                            //num starts
                            s_i[t*sizeSupra + newPacket.source]++;

                            if(DEBUG){
                                sprintf(outputString,"fprintf('NewPacket from %i to %i\\n')",newPacket.source+1,newPacket.destination+1);
                                mexEvalString(outputString);                           
                            }                                
                            //sprintf(outputString,"fprintf('full[0] = %i\\n')",full(&networkQueues[newPacket.source]));
                            //mexEvalString(outputString);                                               

                            //if(!full(&networkQueues[newPacket.source])){
                            if(1){                            
                                //insertedPacketsTemp[t*sizeSupra + newPacket.source] ++;

                                //if(createdPackets < 20){
                                    newPacket.id = createdPackets;
                                    createdPackets = createdPackets + 1;    
                                    //totalPackets[t] = totalPackets[t] + 1;
                                    //BW[newPacket.source] = BW[newPacket.source] + 1;                        
                                    add(&networkQueues[newPacket.source],&newPacket);
                                //}
                            }else{                            
                                sprintf(outputString,"fprintf('Creating: queue of node %i is full. Num elems %i \\n')",
                                        n,numElems(&networkQueues[newPacket.source]));
                                mexEvalString(outputString);  
                                sprintf(outputString,"fprintf('createdPackets = %lu, deliveredPacketsTotal = %lu\\n')",createdPackets,deliveredPacketsTotal);
                                mexEvalString(outputString);                                                                               
                                return;
                            }
                        }else{
                            sprintf(outputString,"fprintf('randIndex = -1, source %i, dest = %i\\n')",n,randDest);
                            mexEvalString(outputString);   
                            return;
                        }
                    }else{
                        if(DEBUG){
                            sprintf(outputString,"fprintf('not possible to go anywhere from = %i, num accessible destination = %i\\n')",n,(int)mxGetN(currAccessibleDest));
                            mexEvalString(outputString);                                  
                        }
                    }
                }                
            }               
        //}
            
        //deliver packets
        for(n=0; n<sizeSupra; n++){
            for(p=0; p<delivRate[n]; p++){                
                
                if(!empty(&networkQueues[n])){
                    
                    newPacket_ptr = consult(&networkQueues[n]);
                    sprintf(outputString,"current node HERE: %d\n", newPacket_ptr->whereAmI);
                    mexEvalString(outputString); 
                    if(newPacket_ptr == NULL && !empty(&networkQueues[n])){
                        sprintf(outputString,"fprintf('possibly trying to consult a virtual element of the queue\\n')");
                        mexEvalString(outputString);     
                        
                        sprintf(outputString,"fprintf('t = %i, numElems = %i, numVirtuals = %i\\n')",
                                t,numElems(&networkQueues[n]),numVirtualElems(&networkQueues[n]));
                        mexEvalString(outputString);     
                        
                        return;
                    } 
                                        
                    if(newPacket_ptr->deliveredOnTime < t){
                                        
                        newPacket_ptr = get(&networkQueues[n]);                                                
                        
                        if(newPacket_ptr == NULL && !empty(&networkQueues[n])){
                            sprintf(outputString,"fprintf('possibly trying to extract a virtual element of the queue\\n')");
                            mexEvalString(outputString);                              
                            return;
                        } 
                        
                        //histogram of packets leaving a node
                        //if(newPacket_ptr->whereAmI != newPacket_ptr->destination){
                        //if(newPacket_ptr->destination != newPacket_ptr->whereAmI){
                        //                if(newPacket_ptr->source != n || newPacket_ptr->maxSteps != -1){                    
                        //if(newPacket_ptr->source != n || newPacket_ptr->maxSteps != -1){                            
                        //if(newPacket_ptr->source == 0){                            
                        //if(newPacket_ptr->previous == 0){                                            
                        //    packetType = newPacket_ptr->source*sizeSupra + newPacket_ptr->destination;                
                        //    distPacketDest[n*sizeSupra*sizeSupra + packetType] ++;                    
                        //}                                                
                        extractedPackets[n] ++;                        
                        
                        extractedPacketsTemp[t*sizeSupra + n] ++;            
                        
                        //if(nextDest == newPacket_ptr->destination || newPacket_ptr->whereAmI == newPacket_ptr->destination){
                        if(newPacket_ptr->whereAmI == newPacket_ptr->destination){
                            
                            //sprintf(outputString,"fprintf('\\nid = %i, delivered = %i, Packet from %i to %i reach destination\\n')",
                            //        newPacket_ptr->id,deliveredPackets,newPacket_ptr->source+1,newPacket_ptr->destination+1);
                            //mexEvalString(outputString); 
                            
                            numPacketOutUnitTime[t] ++;
                            
                            deliveredPacketsTotal = deliveredPacketsTotal + 1;

                            e_i[t*sizeSupra + newPacket_ptr->whereAmI] ++; 
                        }else{
                            
                            //if(n != newPacket_ptr->source){
                                packetsThatContinue[n] ++;
                            //}
                                                        
                            //choose randomly a path
                            // HAUREM DE FER EL MATEIX PERÒ DEPENENT DE L'ORIGEN. ACCEDIR Al MEMBER 'source' DE newPacket_ptr. Llavors PredProb = vector de predProbs [source] i 
                            curCumProb = mxGetCell(PredProb, newPacket_ptr->whereAmI*sizeSupra + newPacket_ptr->destination);// Pred = vector de pred matrix[source]
                            curCumProbPr = mxGetPr(curCumProb);                        
                            curCumPr = mxGetPr(mxGetCell(Pred, newPacket_ptr->whereAmI*sizeSupra + newPacket_ptr->destination));
                            numPred = mxGetM(curCumProb);

                            nextDest = curCumPr[(int)FastGetNextRandomNode_Linear(curCumProbPr, numPred, Uniform(0,1))]-1;                            

                            //sprintf(outputString,"fprintf('NumPred = %i, nextDest = %i\\n')",(int)numPred,nextDest+1);
                            //mexEvalString(outputString);                                                           
                            
                            if(nextDest >= 0 && nextDest < sizeSupra){                                
                              //if(!full(&networkQueues[nextDest])){ 
                              if(1){   
                            
                                //sprintf(outputString,"fprintf('packet from %i to %i (at %i) next = %i\\n')",
                                //        newPacket_ptr->source+1,newPacket_ptr->destination+1,newPacket_ptr->whereAmI+1,nextDest+1);
                                //mexEvalString(outputString);                               
                                  
                                //jumping probabilities
                                //if(n != newPacket_ptr->source){
                                    PExp[newPacket_ptr->whereAmI*sizeSupra + nextDest] ++; 
                                //}                                
                               sigma[t*sizeSupra + nextDest] ++;                                                                    
                                //effective betweenness
                            
                                if(DEBUG){
                                    sprintf(outputString,"fprintf('Packet from from = %i, to = %i, on %i, next dest = %i\\n')",newPacket_ptr->source+1,newPacket_ptr->destination+1,newPacket_ptr->whereAmI+1,nextDest+1);
                                    mexEvalString(outputString);                                                           
                                }
                               
                                //if(newPacket_ptr->whereAmI != newPacket_ptr->source){
                                    if(nextDest != newPacket_ptr->destination){
                                        B_i_ext[t*sizeSupra + nextDest] ++;
                                    }else{
                                        //e_i_ext[t*sizeSupra + nextDest] ++; 
                                    }                                                                  
                                //}
                                if(nextDest != newPacket_ptr->destination){
                                    //B_i[t*sizeSupra + nextDest] ++; 
                                }else{
                                    //e_i[t*sizeSupra + nextDest] ++; 
                                }                                                                                                  
                                  
                                //update deliver time
                                newPacket_ptr->deliveredOnTime = t;
                                newPacket_ptr->previous = newPacket_ptr->whereAmI;                               
                                newPacket_ptr->whereAmI = nextDest;                                
                                newPacket_ptr->maxSteps = 1;
                                
                                if(DEBUG){
                                    sprintf(outputString,"fprintf('New destination = %i\\n')",newPacket_ptr->whereAmI+1);
                                    mexEvalString(outputString);                                   
                                }

                            
                                //insertedPacketsTemp[t*sizeSupra + nextDest] ++;  
                                  
                                //BW[nextDest] = BW[nextDest] + 1;                                
                                                                
                                add(&networkQueues[nextDest],newPacket_ptr);       
                              }else{
                                sprintf(outputString,"fprintf('Routing: queue of node %i is full. Num elems %i \\n')",n,numElems(&networkQueues[nextDest]));
                                mexEvalString(outputString);     
                                sprintf(outputString,"fprintf('createdPackets = %lu, deliveredPacketsTotal = %lu\\n')",createdPackets,deliveredPacketsTotal);
                                mexEvalString(outputString);                                                   
                                return;                              
                              }
                            }else{
                                sprintf(outputString,"fprintf('NumPred = %i, nextDest = %i\\n')",(int)numPred,nextDest);
                                mexEvalString(outputString);   
                                sprintf(outputString,"fprintf('->%i\\n')",nextDest);
                                mexEvalString(outputString);   
                                
                                sprintf(outputString,"fprintf('mexFunction: not possible to jump out of node %i, nextDest = %i\\n')",n,nextDest);
                                mexEvalString(outputString);   
                                sprintf(outputString,"fprintf('Packet from %i to %i\\n')",newPacket.source,newPacket.destination);
                                mexEvalString(outputString);                           
                                return;
                            }
                        }
                    }else{
                        if(DEBUG){
                            sprintf(outputString,"fprintf('\tPacket from %i to %i delivered on time %i\\n')",
                                    newPacket_ptr->source,newPacket_ptr->destination,newPacket_ptr->deliveredOnTime);
                            mexEvalString(outputString);                           
                        }
                    }                    
                }else{
                    //sprintf(outputString,"fprintf('empty queue node %i\\n')",n);
                    //mexEvalString(outputString);                                           
                }
            }                                    
        }   
       
        //numPackets[t] = 0;
        memset ( &numPackets[t*sizeSupra], 0, sizeof(double)*sizeSupra );
        //count packets in the network
        for(n=0; n<sizeSupra; n++){
            numPackets[t*sizeSupra + n] = numElems(&networkQueues[n]);            
        }
        
        if(DEBUG){
            sprintf(outputString,"fprintf('createdPackets = %lu, deliveredPacketsTotal = %lu\\n')",createdPackets,deliveredPacketsTotal);
            mexEvalString(outputString);                   
            sprintf(outputString,"fprintf('\tNum packets t=%i in system %i\\n')",t,(int)numPackets[t]);
            mexEvalString(outputString);               
        }                
    }

    sprintf(outputString,"fprintf('createdPackets = %lu, deliveredPacketsTotal = %lu, genRate[0] = %f\\n')",createdPackets,deliveredPacketsTotal,genRate[0] + genRateProb[0]);
    mexEvalString(outputString);                   

//     //sprintf(outputString,"fprintf('t = %i\\n')",t);
//     //mexEvalString(outputString);                           
//     for(n=0; n<sizeSupra; n++){
//         for(p=0;p<(numElems(&networkQueues[n])-numVirtualElems(&networkQueues[n]));p++){
//             newPacket_ptr = get(&networkQueues[n]);
//             if(newPacket_ptr->whereAmI != newPacket_ptr->previous){
//                 //if(newPacket_ptr->destination != newPacket_ptr->whereAmI){
//                     distPacketDest[n + newPacket_ptr->previous*sizeSupra] ++;                    
//                     //sprintf(outputString,"fprintf('Packets in %i from %i to %i, previous %i\\n')",newPacket_ptr->whereAmI+1,newPacket_ptr->source+1,newPacket_ptr->destination+1,newPacket_ptr->previous+1);
//                     //mexEvalString(outputString);                                               
//                 //}
//             }
//         }
// 
// //         for(p=0; p<sizeSupra; p++){
// //             if(numElems(&networkQueues[n]) > 0){
// //                 //sprintf(outputString,"fprintf('\\tpacketsNode[%i] to %i = %f\\n')",n,p,distPacketDest[n*sizeSupra + p]/numElems(&networkQueues[n]));                
// //                 distPacketDest[n*sizeSupra + p] = distPacketDest[n*sizeSupra + p] / numElems(&networkQueues[n]);
// //             }else{
// //                 //sprintf(outputString,"fprintf('\\tpacketsNode[%i] to %i = %f\\n')",n,p,(double)0);
// //             }
// //             //mexEvalString(outputString);                               
// //         }
//     }
    mxFree(genRate);
    mxFree(delivRate);
    mxFree(networkQueues);
}