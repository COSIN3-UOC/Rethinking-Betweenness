
// to compile: mex -largeArrayDims SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap.cpp -I/usr/local/Cellar/boost/1.80.0/include/
#define TRUE 1
#define FALSE 0
#define DEBUG 0
#include "mex.h"
#define MAX_COLA 2000
#define MAX_NODES 2000
#include "predQueueCircStatic.h"
#include <string.h>
#include <math.h>
#include "matrix.h"
#include <boost/heap/fibonacci_heap.hpp>

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


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    //output matrices
    double tolerance = 0.000000001;
    int numNodes = 0, numLayers = 0;
    int sizeSupra = 0;
    mwSignedIndex dimsP[3];
    
    // get and check inputs
    if (nrhs != 3) mexErrMsgTxt( "Only 3 input arguments allowed." );
    if (nlhs > 5) mexErrMsgTxt( "Only 5 output argument allowed." );
    sizeSupra = mxGetN( prhs[0] );
    if (mxGetM( prhs[0] ) != sizeSupra) mexErrMsgTxt( "Input matrix G needs to be square." );
    if(mxIsSparse(prhs[0])==0) mexErrMsgTxt( "Distance Matrix must be sparse" );

    numNodes = mxGetPr(prhs[1])[0];
    numLayers = mxGetPr(prhs[2])[0];   
    
    if(numLayers != 1){
        mexErrMsgTxt("This code is not prepared to work in multiplex networks" );
    }
    
    
    if((numNodes*numLayers)!=sizeSupra) mexErrMsgTxt( "Incorrect size of the input matrix" );
       
    double *Gpr = mxGetPr(prhs[0]);
    size_t *Gir = mxGetIr(prhs[0]);
    size_t *Gjc = mxGetJc(prhs[0]);   
    //double * GData = mxGetData(prhs[0]);
    double edgeWeight = 0;
    double min_dist = -1;
    
    long int sIndEdges;
    long int eIndEdges;    
    long int nInd;
    int v,s,w,w2,l,temp,t,min_v;
            
    double * C_B = 0;
    double * spEBW_ExcludeInit = 0;
    double * spEBW_Init = 0;
    double * numEndAtNode = 0;
    double * numStartAtNode = 0;
    staticQueue * P = NULL;       
    double * sigma = 0;
    double * sigmaOrig = 0;    
    double * d = 0;
    double * delta = 0;
    double * visited = 0;
    staticQueue S;         
    
    //to store the nodes and the distance to them
    boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> > dk_heap;        
    boost::heap::fibonacci_heap< dk_Node, boost::heap::compare<compare_dk_Node> >::handle_type heapHandles[MAX_NODES];
    
    
    //C_B = (double *)mxMalloc(sizeof(double)*numNodes); 
    plhs[0] = mxCreateDoubleMatrix( sizeSupra, sizeSupra, mxREAL );
    C_B = mxGetPr(plhs[0]);    
    plhs[1] = mxCreateDoubleMatrix( numLayers*numNodes, numLayers*numNodes, mxREAL );
    spEBW_ExcludeInit = mxGetPr(plhs[1]);    
    plhs[2] = mxCreateDoubleMatrix( numLayers*numNodes, numLayers*numNodes, mxREAL );
    spEBW_Init = mxGetPr(plhs[2]);        
    plhs[3] = mxCreateDoubleMatrix( numLayers*numNodes, 1, mxREAL );
    numEndAtNode = mxGetPr(plhs[3]);        
    plhs[4] = mxCreateDoubleMatrix( numLayers*numNodes, 1, mxREAL );
    numStartAtNode = mxGetPr(plhs[4]);                

    visited = (double *)mxMalloc(sizeof(double)*sizeSupra);     
    sigma = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    sigmaOrig = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    d = (double *)mxMalloc(sizeof(double)*sizeSupra); 
    P = (staticQueue *)mxMalloc(sizeof(staticQueue)*sizeSupra);
    delta = (double *)mxMalloc(sizeof(double)*sizeSupra);     
    
    memset ( C_B, 0, sizeof(double)*sizeSupra*sizeSupra );     
    
    for( s = 0 ; s<numNodes ; s++){ 

        sprintf(outputString,"fprintf('Source = %i\\n');",s+1);
        mexEvalString(outputString);
        
        //initializations
        create(&S);                  
        for( w = 0 ; w<sizeSupra ; w++){
            create(&P[w]);       
            sigma[w] = 0;                        
            d[w] = mxGetInf();
            visited[w] = 0;
    
            //store the handle in a table then use the table to access the handle
            heapHandles[w] = dk_heap.push(dk_Node(w,d[w]));             
        }            
        d[s] = 0;
        (*heapHandles[s]).distance = 0; dk_heap.decrease(heapHandles[s]);
        sigma[s] = 1;                
        //end initializations

        //min_v = s;
        //min_dist = 0;       
                
        min_v = dk_heap.top().nodeId;
        min_dist = dk_heap.top().distance;
        dk_heap.pop();        
        
        
        while( min_dist < mxGetInf()){
        //Min_V := S;
        //Min_Dist := 0.0;
        //while Min_Dist /= Plus_Infinity loop

            if(DEBUG){
                sprintf(outputString,"fprintf('min_v %i, min_dist = %f\\n');",min_v+1,min_dist);
                mexEvalString(outputString);                                                             
            }
            
            v = min_v;               
            
            //sprintf(outputString,"fprintf('top = %i\\n');",v+1);
            //mexEvalString(outputString);                                                                     
            
            visited[v] = 1;
            add(&S,v);                                                     
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
                sprintf(outputString,"fprintf('\\tmin_v = %i, min_dist = %f\\n');",min_v+1, min_dist);
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
                    if(d[w] > d[v] + edgeWeight){ // new path to the vertex
                        d[w] = d[v] + edgeWeight;
                        
                        (*heapHandles[w]).distance = d[w]; dk_heap.decrease(heapHandles[w]);                        

                        sigma[w] = sigma[v];
                        create(&P[w]);
                        add(&P[w],v);
                        
                        if(DEBUG){   
                            sprintf(outputString,"fprintf('\\t\\tnew path to vertex %i with dist  = %f\\n');",w+1, d[w]);
                            mexEvalString(outputString);                                                     
                        }
                    }else if( fabs(d[w] - d[v] - edgeWeight) < tolerance){ // the current path is of equal length
                        sigma[w] = sigma[w] + sigma[v];
                        add(&P[w],v);                                                           
                    
                        if(DEBUG){                           
                            sprintf(outputString,"fprintf('\\t\\tpath to equal length to %i , path difference  =%f\\n');",w+1, fabs(d[w] - d[v] - edgeWeight));
                            mexEvalString(outputString);          
                            sprintf(outputString,"fprintf('\\t\\told length = %f, dist to %i = %f and edgeweight = %f\\n');",d[w],v+1,d[v],edgeWeight);
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
        
        //keep an original copy of the sigma
        //memcpy(sigmaOrig,sigma,sizeof(double)*sizeSupra);                
        for( w = 0 ; w<sizeSupra ; w++){
            delta[w] = 0;            
        }          
    
        while(!empty(&S)){
            w = pop(&S);

            if(DEBUG){
                sprintf(outputString,"fprintf('explore vertex w = %i\\n');",w+1);
                mexEvalString(outputString);
            }
                        
            if(DEBUG){
                sprintf(outputString,"fprintf('Distance from (s=%i) to (w=%i) = %f\\n');",s+1,w+1,d[w]);
                mexEvalString(outputString);                                                         
            }
                
            while(!empty(&P[w])){
                v = pop(&P[w]);
            
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\tvertex %i has predecessor %i\\n');",w+1,v+1);
                    mexEvalString(outputString);
                }
                
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\t\\tPred: (v=%i)\\n');",v);
                    mexEvalString(outputString);                                                         
                    sprintf(outputString,"fprintf('\\t\\tsigma[v=%i]=%i/sigma[w=%i]=%i,delta[w=%i]=%f\\n');",v,(int)sigma[v],w,(int)sigma[w],w,delta[w]);
                    mexEvalString(outputString);                           
                }
                
                delta[v] = delta[v] + (sigma[v]/sigma[w])*(1 + delta[w]);                                
                if(v != s){
                    spEBW_ExcludeInit[v*sizeSupra + w] = 
                            spEBW_ExcludeInit[v*sizeSupra + w] + (sigma[v]/sigma[w])*(1 + delta[w]);
                }                
                if(v == s){
                    spEBW_Init[v*sizeSupra + w] = spEBW_Init[v*sizeSupra + w] + (sigma[v]/sigma[w])*(1 + delta[w]);                                                        
                    numStartAtNode[v] = 
                           numStartAtNode[v] + (sigma[v]/sigma[w])*(1 + delta[w]);                    
                }                
                numEndAtNode[w] = numEndAtNode[w] + (sigma[v]/sigma[w]);                                                            
            }
                
            if(w!=s){
                
                C_B[s*sizeSupra + w] = C_B[s*sizeSupra + w] + delta[w];  /// aqui la delta d'ha de normalizar no es poden sumar les betweenness
                                                               // de les differents capes simplement
                
                if(DEBUG){
                    sprintf(outputString,"fprintf('\\t\\tBetweenness[w(%i)]=%f\\n');",w,C_B[s*sizeSupra + w]);
                    mexEvalString(outputString);                                           
                }                
            }
        }
    }
    
    mxFree(sigma);
    mxFree(sigmaOrig);
    mxFree(d);
    mxFree(P);
    mxFree(delta);
    mxFree(visited);
}
