//mex -v -largeArrayDims generate_path_structure.cpp -I"/opt/homebrew/include" -I"/Applications/MATLAB_R2023b.app/extern/include" -L"/opt/homebrew/lib" -L"/Applications/MATLAB_R2023b.app/bin/maca64" -lboost_system -lboost_filesystem -lboost_graph CXXFLAGS="\$CXXFLAGS -std=c++20 -fopenmp" LDFLAGS="\$LDFLAGS -Xpreprocessor -fopenmp"
#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <cstdio>
#include <fstream>
#include <map>
#include <chrono>
#include <numeric>
#include <unordered_map>
#include <future>
#include <mutex>
#include <random>
#include <omp.h>
#include <filesystem>
#include <type_traits>
#include "mex.h"
#include "mat.h"
#include <ctime>
#include <iomanip>
#include <sstream>


#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/copy.hpp>
#include <boost/tokenizer.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/function.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/process.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>



using namespace boost;
using namespace std::chrono;
namespace fs = std::filesystem;

/* ------- LOG and DEBUG ------- */
std::ofstream logfile("app.log", std::ios::app);
std::mutex log_mutex;

std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

void log(const std::string& msg) {
    std::lock_guard<std::mutex> lock(log_mutex);
    logfile << "[" << timestamp() << "] " << msg << std::endl;
}

/* ------------------- GRAPH  DEFINITION -------------------*/

struct VertexProperties {
    int custom_index;
};

typedef boost::adjacency_list<
    boost::vecS,                // OutEdgeList
    boost::vecS,                 // VertexList
    boost::undirectedS,          // UnDirected
    VertexProperties,     // VertexProperties
    boost::property<boost::edge_weight_t, double> // EdgeProperties
> Graph;

// typedef adjacency_list<listS, vecS, undirectedS,
//     no_property, property<edge_weight_t, double> > Graph; 

typedef graph_traits<Graph>::vertex_descriptor Vertex;// descriptors for vertices and edges in the graph, allowing easy access to graph elements.
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef std::pair<Vertex, Vertex> EdgePair;// This represents an edge as a pair of integers (vertex indices).
typedef std::vector<Vertex> Path;
// Define a struct to hold the path property for each vertex
struct VertexProps {
    Path path;
    int id;
};


/* ---------------- TREE DEFINITION ---------------- */


typedef boost::adjacency_list<
    boost::vecS,                // OutEdgeList
    boost::vecS,                 // VertexList
    boost::bidirectionalS,          // Directed
    VertexProps     // VertexProperties
    // boost::no_property // EdgeProperties
> Tree;



typedef graph_traits<Tree>::vertex_descriptor VertexT;// descriptors for vertices and edges in the graph, allowing easy access to graph elements.
typedef graph_traits<Tree>::edge_descriptor EdgeT;
typedef unordered_map<EdgePair, Tree> TreeDict;

/* -------------------- filtered graph type definition -------------------- */

using Filtered = filtered_graph<Graph, boost::function<bool(Edge)>, boost::function<bool(Vertex)> >;

/* ----- HEURISTIC CLASS FOR ASTAR SEARCH ----- */

class Heuristic:
      public boost::astar_heuristic<Filtered, double>
{
public:
    Heuristic(const std::vector<double> distances) : dist_vec(distances) {};

  double operator()(Vertex v) {
    return dist_vec[v];
  }
private:
    const std::vector<double> dist_vec; 

};

/*-- --*/

// exception for astar search termination when goal is found

struct found_goal
{}; 

// visitor that terminates when we find the goal and throws found_goal exception

class astar_goal_visitor : public boost::default_astar_visitor
{
public:
    astar_goal_visitor(Vertex goal) : m_goal(goal) {}
    template < class Filtered > void examine_vertex(Vertex u, Filtered& g)
    {
        if (u == m_goal){
            throw found_goal();}
    }

private:
    Vertex m_goal;
};

/*--------------- FUNCTIONS FOR SPR ALGORITHM ---------------*/


// Concatenates root and spur. Returns tuple with concatenated path and path length.

std::tuple<Path,double> concatenatePaths(Graph G, const Path &path1, const Path &path2) {

    Path concatPath = path1;
    double concatPathLen = 0;

    concatPath.insert(concatPath.end(),
                            path2.begin() + 1, path2.end());

    for (int i = 1; i < concatPath.size(); ++i) {
        Edge e = edge(concatPath[i - 1], concatPath[i], G).first;
        concatPathLen += get(boost::edge_weight, G, e);
    }

    return std::make_tuple(concatPath, concatPathLen);
}

// Finds the spur path through an ASTAR SEARCH from the spur node to the target according to Yen's alg.

std::tuple<Path, double> FindShortestPathYen(Graph &G, int spurNode, int target,
                                            const std::vector<EdgePair> &devLinks,
                                            const Path &rootPath, std::vector<double>& costsPaths){

    // Afiltered "copy" f with the removed nodes from the root path and deviation links:
    std::set<Vertex> removed_set;
    std::set<Edge> removed_edge_set;
    Filtered f(G, [&](Edge e){ return removed_edge_set.end() == removed_edge_set.find(e); },
            [&](Vertex v){ return removed_set.end() == removed_set.find(v); });

    Path spurPath;


    // Remove deviation links from the graph
    for (const auto &link : devLinks) {
        // remove_edge(link.first, link.second, GOrig);
        removed_edge_set.insert(boost::edge(link.first, link.second,G).first);
    }

    // Remove nodes from the root path
    for (int i = 0; i < rootPath.size() - 1; ++i) {
        Vertex u = rootPath[i];
        removed_set.insert(u);
        // clear_vertex(u, GOrig);
    }

    int spurDeg = degree(spurNode, f);
            if (spurDeg == 0){
                return std::make_tuple(spurPath, std::numeric_limits<double>::max());
            };


    int targDeg = degree(target, f);
    if (targDeg == 0){
        return std::make_tuple(spurPath, std::numeric_limits<double>::max());
    };

    // Dijkstra's/A* algorithm to find the spur path

    std::size_t numVertices = std::distance(boost::vertices(G).first, boost::vertices(G).second);

    Path predecessors(numVertices);
    std::vector<double> distances(numVertices);

    Heuristic heuristic(costsPaths);
    astar_goal_visitor vis(target);

    try
    {
    boost::astar_search(f, spurNode, heuristic,
                        predecessor_map(make_iterator_property_map(predecessors.begin(), get(vertex_index, f)))
                        .distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, f)))
                        .visitor(vis));
    } 
    catch (found_goal fg){

        // Reconstruct the spur path

    int currentNode = target;

    while (currentNode != spurNode) {

        spurPath.emplace_back(currentNode);
        if (predecessors[currentNode] == currentNode) {
            spurPath.clear();
            return std::make_tuple(spurPath, std::numeric_limits<double>::max());
        }
        currentNode = predecessors[currentNode];
    }
    spurPath.emplace_back(spurNode);
    std::reverse(spurPath.begin(), spurPath.end());
    }

    // clearing removed sets
    removed_set.clear();
    removed_edge_set.clear();

    return std::make_tuple(spurPath, distances[target]);

}

//Gets the candidate spur path from the initial shortest path tree (Dijkstra from target) from the VF algorithm

Path GetCandidatePathVF(const Graph& G, std::vector<EdgePair> dev_links, 
                        int dev_node, std::vector<double>& costs_paths, 
                        Path& possible_paths,
                        int targ) {

    Path cand_path;
    double cost = std::numeric_limits<int>::max();
    Path cand_path_temp;
    double candidate_cost;

// For each out edge from the spur node dfferent than the deviation link
// get the spur path stored in the initial shortest path tree and get the shortest
// spur path 

    for (auto succ : make_iterator_range(out_edges(dev_node, G))) {
        EdgePair e = EdgePair(source(succ, G), target(succ, G));
        if (std::find(dev_links.begin(), dev_links.end(), e) == dev_links.end()) {

            if (source(succ, G) == dev_node){
                candidate_cost = get(boost::edge_weight, G, succ) + costs_paths[target(succ, G)];
                Vertex currentNode = target(succ, G);

                while (currentNode != targ) {
                    cand_path_temp.emplace_back(currentNode);
                    currentNode = possible_paths[currentNode];
                }
                cand_path_temp.emplace_back(targ);
            }
            else{
                candidate_cost = get(boost::edge_weight, G, succ) + costs_paths[source(succ, G)];
                
                Vertex currentNode = source(succ, G);

                while (currentNode != targ) {
                    cand_path_temp.emplace_back(currentNode);
                    currentNode = possible_paths[currentNode];
                }
                cand_path_temp.emplace_back(targ);
            }
            if (candidate_cost < cost) {
                cand_path = cand_path_temp;
                cost = candidate_cost;
            }

            cand_path_temp.clear();

        }
    }

    return cand_path;
    }

    // Function to retrieve the candidate path from the graph.
    //Iterates through the stored path tree and looks for a valid spur path (no intersection with root path).
    // If no valid spur path is found it returns the last path stored in the tree as candidate

    std::tuple<Path, Path, VertexT> RetrieveCandidatePath(Tree& candTree, const std::vector<Vertex>& N, Vertex devNode) {
    Path spurPath;
    Path candPath;
    VertexT lastNode = graph_traits<Tree>::null_vertex();

    std::queue<Vertex> nodeQueue;
    nodeQueue.push(vertex(0, candTree)); 

    while (!nodeQueue.empty()) {

        VertexT currentNode = nodeQueue.front();
        nodeQueue.pop();

        candPath = candTree[currentNode].path;

        // Check if the candidate path is a valid spur path

        std::set<Vertex> removedSet(N.begin(), N.end());
        removedSet.insert(devNode);

        // Check intersection with root to see if it is a valid spur path

        bool isSpurPath = std::all_of(candPath.begin(), candPath.end(), [&](Vertex node) {
            return removedSet.find(node) == removedSet.end();
        });
        if (isSpurPath) {
            lastNode = currentNode;

            spurPath = candPath;
            break;
        }

        // Add children to the queue
        auto childrenRange = out_edges(currentNode, candTree);

        for (auto it = childrenRange.first; it != childrenRange.second; ++it) {
            if (std::find(N.begin(), N.end(), candTree[target(*it, candTree)].id) != N.end()){

                nodeQueue.push(target(*it, candTree));
            }
        }
        lastNode = currentNode;
    }
    return std::make_tuple(spurPath, candPath, lastNode);
}


// Function to calculate candidate paths. Grows the candidate path tree and finds the spur path.

Path CalculateCandidatePaths(Graph& G, Tree& candTree, const Path& N,
                            Path& candPath, VertexT lastLeafNode,
                            const std::pair<Vertex, Vertex>& devLink, Vertex devNode, Vertex source,
                            Vertex target, std::vector<double>& costsPaths) {

    // Step 1: Initialization
    VertexT next_node;
    VertexT parent_node = lastLeafNode;
    VertexT current_node = lastLeafNode;
    boost::graph_traits<Tree>::vertex_iterator vi, vi_end;
    // filtered graph according to the candidate path tree.
    std::set<Vertex> removed_set;
    std::set<Edge> removed_edge_set;
    Filtered f(G, [&](Edge e){ return removed_edge_set.end() == removed_edge_set.find(e); },
        [&](Vertex v){ return removed_set.end() == removed_set.find(v); });

    if (num_vertices(candTree) > 1){
    
        // from the last leaf of the tree, remove all previous vertices from the graph
        while (current_node != graph_traits<Tree>::null_vertex()) {
            //remove node from filtered graph

            if (candTree[current_node].id != -1){
                removed_set.insert(candTree[current_node].id);
            }

            auto inEdges = in_edges(current_node, candTree);

            // if no in edges
            if (inEdges.first == inEdges.second){
                next_node = 0;
                break;
            }

            for (auto ei = inEdges.first; ei != inEdges.second; ++ei) {
                next_node = boost::source(*ei, candTree);
            }

            if (next_node != current_node) {
                current_node = next_node;
            } 
            else {
                break;
            }
        }
    }
    // removing devition link
    if (edge(devLink.first, devLink.second, f).second) {
    removed_edge_set.insert(edge(devLink.first, devLink.second, G).first);
    }

    // Step 2: Candidate path calculation

    int rm_node = -1;

    std::set<Vertex> removed_node_set(N.begin(), N.end());
    removed_node_set.insert(devNode);

    // While the intersection between the candidate spur and root paths is not 0 keep expanding the tree

    while (!candPath.empty() && !std::all_of(candPath.begin(), candPath.end(),
                                    [&](Vertex node) {return removed_node_set.find(node) == removed_node_set.end(); })) {

    // Select node nearest to deviation node in candPath intersection N
    Path intersec;
    Path sorted_cand = candPath;
    Path sorted_N = N;

    std::sort(sorted_cand.begin(), sorted_cand.end());
    std::sort(sorted_N.begin(), sorted_N.end());

    std::set_intersection(sorted_cand.begin(), sorted_cand.end(), sorted_N.begin(), sorted_N.end(), std::back_inserter(intersec));

    std::map<Vertex, Vertex> dists;
    if (!intersec.empty()){   
        for (Vertex node : intersec) {

            Vertex n1 = devNode;
            double total_dist = 0;

            for (int i = N.size() - 1; i >= 0; --i) {
                Vertex n2 = N[i];
                if (edge(n1, n2, f).second) {
                    total_dist += get(edge_weight, f, edge(n1, n2, f).first);
                } else if (edge(n2, n1, f).second) {
                    total_dist += get(edge_weight, f, edge(n2, n1, f).first);
                }
                n1 = n2;
            }
            dists[node] = total_dist;
        }

        auto min_dist_node = std::min_element(dists.begin(), dists.end(),
                                        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        rm_node = min_dist_node->first;
        // Remove rm_node
        removed_set.insert(rm_node);
    }

    // Use dijkstra/astar to find shortest path
    Path p(num_vertices(G));
    std::vector<double> d(num_vertices(G));
    Vertex s = devNode;

    Heuristic heuristic(costsPaths);
    astar_goal_visitor vis(target);

    if (degree(s, f) == 0 or degree(target, f) == 0){
        Path empt;
        return empt;
    }
    try
    {

        boost::astar_search(f, s, heuristic,
                        predecessor_map(make_iterator_property_map(p.begin(), get(vertex_index, f)))
                        .distance_map(make_iterator_property_map(d.begin(), get(vertex_index, f)))
                        .visitor(vis));
        
        // if the target node is not reachable from the deviation node return an empty path
        Path empt;
        return empt;
    } 

    catch (found_goal fg){

        // Reconstruct the spur path

        candPath.clear();
        Vertex prev_v = graph_traits<Graph>::null_vertex();

        for (Vertex v = target; v != s; v = p[v]) {
            candPath.emplace_back(static_cast<Vertex>(v));
            if (prev_v == v){
                Path empt;
                return empt;
            }
            prev_v = v;
        }

    reverse(candPath.begin(), candPath.end());   
    }


    // Update the candidate tree with the new path

    if (rm_node != -1) {

        VertexT new_node = add_vertex(candTree);
        candTree[new_node].id = rm_node;
        candTree[new_node].path = candPath;
        add_edge(parent_node, new_node, candTree);
        parent_node = new_node;


    } else {
        candTree[parent_node].id = devNode;
        candTree[parent_node].path = candPath;
    }
    }
    removed_set.clear();
    removed_edge_set.clear();
    return candPath;
}

// uses all the above functions to find the spur path from a deviation node given a set of dev. links

std::pair<Path, int> FindSpurPath_SPR(Graph& G, const std::vector<EdgePair>& devLinks, const Path& N, Vertex devNode,
                std::vector<double>& costsPaths, Path& possiblePaths,
                TreeDict& treeDict, const Path& rootPath, Vertex source, Vertex target) {

    int reused = 0;
    Path spurPath;
    Path candPath;
    Tree candTree;

    if (devLinks.size() == 1) { // SPR procedure only if there is 1 deviation link
        EdgePair devLink = devLinks[0];

        try {
            candTree = treeDict.at(devLink);
            // If the candidate tree is already stored
        } catch (const std::out_of_range&) {
            // If the candidate tree is not stored, retrieve it
            candPath = GetCandidatePathVF(G, devLinks, devNode, costsPaths, possiblePaths, target);

            //create tree
            Vertex u = boost::add_vertex(candTree);
            candTree[u].path = candPath;
            candTree[u].id = -1;

            //add tree to the dict
            treeDict[devLink] = candTree;
        }

        //determining spur (if there is) and candidate paths

        Path iniCandPath;
        Vertex lastLeafNode;

        std::tie(spurPath, iniCandPath, lastLeafNode) = RetrieveCandidatePath(candTree, N, devNode);   


        //with the VF and tree we could reuse a previous path
        if (!spurPath.empty()) {
            spurPath.insert(spurPath.begin(), devNode);
            reused++;
            
            // printf("reused paths: %d\n", reused);
            return std::make_pair(spurPath, reused);
        } 
        
        else { //no reuse available

            spurPath = CalculateCandidatePaths(G, candTree, N, iniCandPath, lastLeafNode, devLink, 
                                                devNode, source, target, costsPaths);
            if (spurPath.empty()) {
                // No spur paths exist given this root and target
                return std::make_pair(spurPath, reused);
            } else {
                spurPath.insert(spurPath.begin(), devNode);
                treeDict[devLink] = candTree;

                //  printf("dev Node %d", devNode);
                return std::make_pair(spurPath, reused);
            }
        }

    } else {
        candPath = GetCandidatePathVF(G, devLinks, devNode, costsPaths, possiblePaths, target);

    std::set<Vertex> intersection;
    intersection.insert(N.begin(), N.end());
    intersection.insert(devNode);

    if (std::all_of(candPath.begin(), candPath.end(), [&](Vertex node) { return intersection.find(node) == intersection.end(); })) {
        candPath.insert(candPath.begin(),devNode);
        spurPath = candPath;
        reused++;
        return std::make_pair(spurPath, reused);
        // spur_path_length = cand_cost;
    }else {
            std::tuple<Path, double> spurPath_t = FindShortestPathYen(G, devNode, target, devLinks, rootPath, costsPaths);
            return std::make_pair(std::get<0>(spurPath_t), reused);

            }
    }
}
// Finds deviation links given previous paths and root path
std::vector<EdgePair> findDeviationLinks(Graph G, std::vector<Path> paths, Path rootPath, int sizeOfRoot){

    std::vector<EdgePair> devLinks;
    for (Path ipath : paths) {
        if (ipath.size() > sizeOfRoot + 1){
            Path prev_path(ipath.begin(), ipath.begin() + sizeOfRoot + 1);

            if (rootPath == prev_path){
                //finding deviation links
                Vertex u = ipath[sizeOfRoot];
                Vertex v = ipath[sizeOfRoot + 1];

                //  may not be needed ---------------------------------
                //check if dev. link is in graph
                bool exists = false;
                std::pair<Edge, bool> edgeInfo = edge(u, v, G);
                exists = edgeInfo.second;
                // ----------------------------------------------------
                EdgePair e_pair = std::make_pair(u, v);
                if (exists && std::find(devLinks.begin(), devLinks.end(), e_pair) == devLinks.end()) {
                    devLinks.emplace_back(e_pair);
                } 
            }
        }
    }
    return devLinks;
}

/*

************************ SPR MAIN FUNCTION **************************

*/

std::tuple<std::vector<Path>, std::vector<double>> TolShortestPaths_SPR(Graph G, Vertex source, Vertex target, Path& ini_predecessors, 
                                                                        std::vector<double>& distances, const double tol, bool debug=false) {
    struct CompareTuples {
        bool operator()(const std::tuple<double, Path>& a, const std::tuple<double, Path>& b) const {
        return std::get<0>(a) > std::get<0>(b); // Greater comparison for min-heap based on the first element of the tuple (weight)
        }
    };

    std::vector<Path> real_paths;
    std::vector<Path> paths;
    std::vector<double> pathLens;
    TreeDict treeMap;
    std::vector<std::tuple<double, Path>> B;
    std::make_heap(B.begin(), B.end(), CompareTuples());
    paths.clear();

    // INITIAL DIJKSTRA -----

    if (ini_predecessors.empty()){
        
        ini_predecessors.resize(num_vertices(G), -1);  // Initialize path with -1 for all vertices
        distances.resize(num_vertices(G), 0); // Store distances for Dijkstra's algorithm

        dijkstra_shortest_paths(G, target,
                                predecessor_map(make_iterator_property_map(ini_predecessors.begin(), get(vertex_index, G)))
                                .distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, G))));

    }

    // Extract the path to the target node
    Path pathToTarget;

    for (Vertex v = source; v != target; v = ini_predecessors[v]) {
        
        pathToTarget.emplace_back(v);

    }

    pathToTarget.emplace_back(target);
    double dist_SP = distances[source]; //distance of the shortest path

    // -----

    real_paths.emplace_back(pathToTarget);
    pathLens.emplace_back(dist_SP);    

    // Compute another dijkstra to know the distances to the source node
    Path preds;
    std::vector<double> dists_to_s;

    preds.resize(num_vertices(G), -1);  // Initialize path with -1 for all vertices
    dists_to_s.resize(num_vertices(G), 0); // Store distances for Dijkstra's algorithm

    dijkstra_shortest_paths(G, source,
                            predecessor_map(make_iterator_property_map(preds.begin(), get(vertex_index, G)))
                            .distance_map(make_iterator_property_map(dists_to_s.begin(), get(vertex_index, G))));


    std::set<Vertex> outranged_n;
    std::set<Edge> outranged_e;


    Filtered G_red(G, [&](Edge e){ return outranged_e.end() == outranged_e.find(e); },
            [&](Vertex v){ return outranged_n.end() == outranged_n.find(v); });

    // Iterate through the vector using a range-based for loop
    double epsilon = pow(10,-10);
    for (int i = 0; i < dists_to_s.size(); ++i) {
        double tot_ds = distances[i] + dists_to_s[i];
        if (tot_ds > dist_SP*(1+tol)) {

            //check if the difference in values is actually significant
            if (std::abs(tot_ds - dist_SP*(1+tol)) > epsilon){

                outranged_n.insert(i);
            }
            
        }
    }       

    Graph reduced_G;
    copy_graph(G_red, reduced_G);

    // printf("reduced total number of nodes from %zu to %zu\n", num_vertices(G), num_vertices(reduced_G));

    // Access the vertex index property map
    auto vertexIndexMap = get(boost::vertex_index, reduced_G);

    // map of old ID to new ID
    std::unordered_map<std::size_t, Vertex> inverseMap;
    boost::graph_traits<Graph>::vertex_iterator vertexIt, vertexEnd;

    for (boost::tie(vertexIt, vertexEnd) = boost::vertices(reduced_G); vertexIt != vertexEnd; ++vertexIt) {
        std::size_t vertexID = boost::get(vertexIndexMap, *vertexIt);
        std::size_t nodeID = reduced_G[*vertexIt].custom_index;
        inverseMap[nodeID] = vertexID;
    }

    // New source and target IDs

    Vertex new_s = inverseMap[source];
    Vertex new_t = inverseMap[target];
    typedef boost::graph_traits<Filtered>::edge_iterator edge_iterator;
    
    std::pair<edge_iterator, edge_iterator> edges = boost::edges(G_red);
    
    //Do the initial dijkstra with the new node ids

    Path ini_preds;
    std::vector<double> dists;
    ini_preds.resize(num_vertices(reduced_G), -1);  // Initialize path with -1 for all vertices
    dists.resize(num_vertices(reduced_G), 0); // Store distances for Dijkstra's algorithm

    dijkstra_shortest_paths(reduced_G, new_t,
                            predecessor_map(make_iterator_property_map(ini_preds.begin(), get(boost::vertex_index, reduced_G)))
                            .distance_map(make_iterator_property_map(dists.begin(), get(boost::vertex_index, reduced_G))));

    

    // Extract the path to the new target node id ------
    Path new_SP;

    for (Vertex v = new_s; v != new_t; v = ini_preds[v]) {        
        new_SP.emplace_back(v);
    }

    new_SP.emplace_back(new_t);
    // ------

    paths.emplace_back(new_SP);

    double l = dist_SP;
    // int tot_spur = 0;
    // int n_reused = 0;
    int itercount = 0;
    while (l <= dist_SP * (1 + tol)) {

        if (itercount%5000 == 0 and debug and itercount > 1){

            // printf("Number of paths: %zu, " 
            // "Path len and max len: %lf %lf, " 
            // "(s, t): (%d, %d)\n", 
            // paths.size(), l, dist_SP * (1 + tol), source, target);

            char out_log[100];
            // Writing successors
            std::sprintf(out_log, "Number of paths: %zu, " 
            "Path len and max len: %lf %lf, " 
            "(s, t): (%d, %d)\n", 
            paths.size(), l, dist_SP * (1 + tol), source, target);
            log(out_log);
        }
        for (int j = 0; j < paths.back().size() - 1; ++j) {
            Vertex spurNode = paths.back()[j];

            Path rootPath(paths.back().begin(), paths.back().begin() + j + 1);//takes the root path from the last path using iterators

            //find deviation links
            std::vector<EdgePair> devLinks = findDeviationLinks(reduced_G, paths, rootPath, j);

            // if the degree of the spur node is the same as the amount of deviation 
            // links continue with the next spur node
            int spurDeg = degree(spurNode, reduced_G);

            if (spurNode == new_s){

                if (spurDeg == devLinks.size()){

                    continue;
                };
                
            } else {

                if (devLinks.size() >= (spurDeg-1)){

                    continue;
                };
            }

            //nodes that have to be removed from the graph (root path)

            Path N(paths.back().begin(), paths.back().begin() + j);
            // tot_spur++;

            //getting spur path
            auto spurPath_tp = FindSpurPath_SPR(reduced_G, devLinks, N, spurNode, dists, 
                                            ini_preds, treeMap, rootPath, new_s, new_t);  
            

            Path spurPath = std::get<0>(spurPath_tp);
            // n_reused += std::get<1>(spurPath_tp);
            
            if (spurPath.size() == 0){

                std::tuple<Path, double> spurPath_t = FindShortestPathYen(reduced_G, spurNode, new_t, devLinks, rootPath, dists);
                spurPath = std::get<0>(spurPath_t); 

            }
            
            if (!spurPath.empty() && spurPath.back() == new_t) 
                { 
                
                //concatenating and retrieving total path and total path length
                std::tuple totalPath_t = concatenatePaths(reduced_G, rootPath, spurPath); 

                Path totalPath = std::get<0>(totalPath_t);
                double totalPathLength = std::get<1>(totalPath_t);

                //add path only if its length is relevant otherwise don't even bother
                if (totalPathLength <= dist_SP * (1 + tol)){

                    //finding if total path is in the heap and adding it accordingly.
                    auto it = std::find_if(B.begin(), B.end(), [&totalPath](const auto& element) {
                    return std::get<1>(element) == totalPath;
                    });
                    spurPath.clear();

                    if (it == B.end() || B.size() == 0) {
                        B.emplace_back(std::make_pair(totalPathLength, totalPath));
                        std::push_heap(B.begin(), B.end(), CompareTuples());

                }
                }
            }
        }


        if (!B.empty()) {
            auto lp = B.front();
            std::pop_heap(B.begin(), B.end(), CompareTuples());
            B.pop_back();
            l = std::get<0>(lp);
            auto p = std::get<1>(lp);

            if (l <= dist_SP * (1 + tol)) {
                // printf("closeness to Tolerance = %lf\n", (l-dist_SP)/(dist_SP * tol));
                paths.emplace_back(p);

                Path real_p;
                for (int i = 0; i < p.size(); ++i) {
                    real_p.emplace_back(reduced_G[p[i]].custom_index);
                }
                real_paths.emplace_back(real_p);
                pathLens.emplace_back(l);

            }
        } else {
            break;
        }
        itercount++;
    }
    // printf("reused ratio %lf\n", n_reused/double(tot_spur));
    return std::make_tuple(real_paths, pathLens);
}


// Conversion function from std::vector<unsigned long> to mxArray*

template<typename T>
mxArray* convertVectorToMxArray(const std::vector<T>& vec) {
    mwSize dims[] = {vec.size()};
    mxClassID classID;

    // Determine the appropriate MATLAB class ID based on the C++ data type
    if constexpr (std::is_same<T, float>::value) {
        classID = mxSINGLE_CLASS;
    } else if constexpr (std::is_same<T, double>::value) {
        classID = mxDOUBLE_CLASS;
    } else if constexpr (std::is_same<T, int32_t>::value) {
        classID = mxINT32_CLASS;
    } else if constexpr (std::is_same<T, Vertex>::value) {
        classID = mxINT32_CLASS;
    } else {
        // Add more types as needed
        // mexErrMsgIdAndTxt("convertVectorToMxArray:unsupportedType", "Unsupported data type for conversion.");
        return nullptr; // For safety, in case mexErrMsgIdAndTxt returns
    }

    mxArray* mxVec = mxCreateNumericArray(1, dims, classID, mxREAL);
    if (!mxVec) {
        // mexErrMsgIdAndTxt("convertVectorToMxArray:creationFailed", "Failed to create MATLAB array.");
        return nullptr;
    }

    void* mxData = mxGetData(mxVec);
    if (!mxData) {
        // mexErrMsgIdAndTxt("convertVectorToMxArray:accessFailed", "Failed to access data of MATLAB array.");
        mxDestroyArray(mxVec); // Cleanup before returning
        return nullptr;
    }

    // Copy data from the input vector to the mxArray
    std::copy(vec.begin(), vec.end(), static_cast<T*>(mxData));

    return mxVec;
}

// write the cell matrix
void writeCellMatrixToFile(const char* filename, const char* variable_name, mxArray* cellMatrix) {
    MATFile *pmatFile;
    
    // Open file for writing
    pmatFile = matOpen(filename, "w");
    // if (pmatFile == NULL) {
    //     mexErrMsgIdAndTxt("MyProg:FileWriteError", "Error opening file %s for writing.\n", filename);
    //     return;
    // }
    // Write the cell matrix to the file
    if (matPutVariable(pmatFile, variable_name, cellMatrix) != 0) {
        // mexErrMsgIdAndTxt("MyProg:FileWriteError", "Error writing cell matrix to file %s.\n", filename);
        matClose(pmatFile);
        return;
    }

    // Close the file
    if (matClose(pmatFile) != 0) {
        // mexErrMsgIdAndTxt("MyProg:FileCloseError", "Error closing file %s.\n", filename);
        return;
    }
}

/*------------------------------------TOLERANCE BETWEENNESS------------------------------------*/
Graph adj_mat_to_graph(const mxArray *A){

    // Ensure the input is a 2D array
    // if(mxGetNumberOfDimensions(A) != 2) {
    //     mexErrMsgIdAndTxt("MATLAB:conversion", "Input must be a 2D array.");
    // }

    // Get the dimensions of the mxArray
    size_t rows = mxGetM(A);
    size_t cols = mxGetN(A);

    // Get a pointer to the data in mxArray
    double* data = mxGetPr(A);

    // Initialize the 2D std::vector
    std::vector<std::vector<float>> Avec(cols, std::vector<float>(rows));

    Graph G;
    // Copy data from the mxArray to the std::vector
    // Remembering that MATLAB uses column-major order
    for (size_t col = 0; col < cols; ++col) {
        for (size_t row = 0; row < rows; ++row) {

            Edge e = add_edge(row, col, G).first;
            put(edge_weight, G, e, static_cast<float>(data[col * rows + row]));
        }
    }

return G;
}

//template<typename doub_vec>
// void compute_tolSP_paths(const mxArray *A, doub_vec tols) {

/**
 * @brief Compute all paths in a graph that are within a given tolerance of the shortest path
 *        for every source-target pair and save the results to a MATLAB cell array file.
 * 
 * This function performs the following steps:
 * 1. For each source-target pair (s < t), computes all paths whose lengths are within
 *    a tolerance `tol` of the shortest path.
 * 3. Stores each set of paths in a MATLAB MxN matrix, where each row corresponds to a path,
 *    padded to the maximum path length (number of nodes in the network), and the last column stores the path length.
 * 4. Saves the full NxN cell array of paths to a `.mat` file.
 * 
 * @param G Reference to the input graph object (Boost Graph Library type).
 * @param tol Tolerance factor for path length relative to the shortest path. 
 *            Paths longer than (1+tol) * shortest_path_length are discarded.
 * @param graph_number An integer identifier used to generate the output file name.
 * 
 * Notes:
 * - Uses OpenMP to parallelize the computation over source-target pairs.
 * - MATLAB indices are 1-based; node indices in the matrix are incremented by 1.
 * - Each cell of the output MATLAB cell array contains a matrix of type double:
 *      * Rows: individual paths from source to target
 *      * Columns: node indices (padded if necessary) + last column as path length
 * - The output file is stored under `./ring_road_p_structure_vx10/ring_road_<graph_number>/path_ls_mat_<tol>.mat`.
 */

void compute_tolSP_paths(Graph &G, double tol, int graph_number, bool debug){

    int nnodes = num_vertices(G);
    
    std::unordered_map<Vertex, std::tuple<Path, std::vector<double>>> dijkstra_dict;
    
    // initial dijkstra
    Path ini_predecessors(num_vertices(G));
    std::vector<double> ini_dists(num_vertices(G));
    for (size_t i = 0; i < nnodes; ++i) {
        dijkstra_shortest_paths(G, i,
                        predecessor_map(make_iterator_property_map(ini_predecessors.begin(), get(vertex_index, G)))
                        .distance_map(make_iterator_property_map(ini_dists.begin(), get(vertex_index, G))));

        dijkstra_dict[i] = std::make_tuple(ini_predecessors, ini_dists);
    }

    std::tuple<Path, std::vector<double>> ini_path_t;
    int it = 0;

    std::unordered_map<int, std::tuple<int, int>> loop_map;

    for (int i = 0; i < nnodes - 1; ++i) {
        for (int j = i + 1; j < nnodes; ++j) {

            loop_map[it] = std::make_tuple(i, j);
            it++;
        }
    }

    // Create a cell array to store the results
    mxArray *cellPathMat = mxCreateCellMatrix(nnodes, nnodes);

    // beginning parallel loop
    it = 0;
    // Preallocate a vector for each thread
    std::vector<mxArray*> thread_results(nnodes*nnodes, nullptr);

    double last = omp_get_wtime();

    #pragma omp parallel private(ini_path_t) shared(it, tol, thread_results, debug, last)
    {
        double last = omp_get_wtime();
        #pragma omp for schedule(dynamic,1) 
        for (int k = 0; nnodes*(nnodes-1)/2 > k; k++) {

            auto [s, t] = loop_map[k];
            ini_path_t = dijkstra_dict.at(t);

            // compute the QSPs

            auto pair_t = TolShortestPaths_SPR(G, s, t, std::get<0>(ini_path_t), std::get<1>(ini_path_t), tol, debug);
    
            std::vector<double> path_lens = std::get<1>(pair_t);

            auto paths = std::get<0>(pair_t);
            // Suppose paths is std::vector<Path>, Path = std::vector<int>
            size_t num_paths = paths.size();
            size_t L = nnodes+1; // could also vary per path!

            // Create an num_paths-by-L MATLAB matrix of type double
            mxArray *mat = mxCreateDoubleMatrix(num_paths, L, mxREAL);
            double *data = mxGetDoubles(mat);

            // Fill it row by row
            for (size_t i = 0; i < num_paths; i++) {
                for (size_t j = 0; j < paths[i].size(); j++) {
                    data[i + j*num_paths] = paths[i][j]; 
                }

                // add the path length at the end of the row
                data[i + num_paths*(L-1)] = path_lens[i]; 
            }

            int idx = s + t*nnodes; // column-major order (row=s, col=t)
            thread_results[idx] = mat;
            #pragma omp atomic
            it++;

            // print occasionally
            // if (it % 100 == 0) {
            //     #pragma omp critical
            //     {
            //         printf("Progress: %.4f\n", 2*it/double(nnodes*(nnodes-1)));
            //     }
            // }  

            
            if (omp_get_wtime() - last > 120) { // once per 2 mins
                #pragma omp critical
                {
                printf("Graph %d, Progress: %.4f\n", graph_number, 2*it / double(nnodes*(nnodes-1)));
                char out_log[100];
                // Writing successors
                sprintf(out_log, "Progress: %.4f\n Using %d threads\n",
                2*it / double(nnodes*(nnodes-1)), omp_get_num_threads());
                log(out_log);
                last = omp_get_wtime();
                }
            }
        }
    }
        // After parallel loop (single-threaded)
        for(int i = 0; i < nnodes*nnodes; ++i) {
            if(thread_results[i] != nullptr)
                mxSetCell(cellPathMat, i, thread_results[i]);
        }
        //write cell array
        char fname[100];
        char fname2[100];
        char path_ls_nm[30];
        // Writing successors
        std::sprintf(fname, "./star_road_p_structure_vx2/ring_road_%d", graph_number);
        boost::filesystem::create_directories(fname);
            
        std::sprintf(fname2, "/path_ls_mat_%d.mat", static_cast<int>(tol*100));

        std::string file_dir = std::string()+fname+fname2;
        const char* cString = file_dir.c_str();

        std::sprintf(path_ls_nm, "path_list");//%d", orig);
        writeCellMatrixToFile(cString, path_ls_nm, cellPathMat);
        
}

int extractFirstNumber(const std::string& filename) {
    std::string number;
    for (char ch : filename) {
        if (std::isdigit(ch)) {
            number += ch;
        } else if (!number.empty()) {
            // We've reached the end of a number
            break;
        }
    }

    return number.empty() ? -1 : std::stoi(number);  // Return -1 if no number found
}

int main() {

    /*
    
    READING CSV FILE: EDGE[0], EDGE[1], WEIGHT

    */

   fs::path graph_folder = "./star_roads/vx2/graph_variations";

   // List the directories in the graph_folder
    std::vector<fs::directory_entry> entries;
    for (const auto& e : fs::directory_iterator(graph_folder)) entries.push_back(e);
    // std::sort(entries.begin(), entries.end(), [](auto& a, auto& b){ return a.path() < b.path(); });

    const int start = 0;  // skip the first 4
    omp_set_max_active_levels(2);   // or export OMP_MAX_ACTIVE_LEVELS=2
    omp_set_nested(1);              // or export OMP_NESTED=TRUE
    int tot_size = (int)entries.size();
    int completed = 0;
    // 2) Parallelize with teams + (optional) inner parallel
    #pragma omp parallel num_threads(3)
    {
        // distribute files to 2 outer threads
        #pragma omp for schedule(static, (entries.size()+ 2)/3)
        for (int file_count = start; file_count < (int)entries.size(); ++file_count) {
            const auto& entry = entries[file_count];
            if (entry.path().extension() == ".csv") {
                // Read the graph from the CSV file
                fs::path graph_fname =  entry.path();
                std::string filename = entry.path().filename().string();  // e.g., "file_123.csv"
                int graph_number = extractFirstNumber(filename);

                // log print
                printf("Processing file: %d/%d graph number %d\n", file_count, tot_size, graph_number);
                char out_log[100];
                std::sprintf(out_log, "Processing file: %d/%d graph number %d\n", file_count, tot_size, graph_number);
                log(out_log);

                Graph read_g(0);
                std::ifstream file(graph_fname);
                std::string line;

                // Header ignored (1st line)
                std::getline(file, line);

                while (std::getline(file, line)) {
                    boost::tokenizer<boost::escaped_list_separator<char>> tokens(line);
                    auto tokenIterator = tokens.begin();
                    Vertex vertex1 = std::stoi(*tokenIterator++);
                    Vertex vertex2 = std::stoi(*tokenIterator++);
                    double weight = std::stod(*tokenIterator);

                    // Add edge to the graph with the given weight
                    Edge e = add_edge(vertex1, vertex2, read_g).first;
                    put(edge_weight, read_g, e, weight);
                }

                // obtain a property map for the custom_index property
                property_map<Graph,int VertexProperties::*>::type index = get(&VertexProperties::custom_index, read_g);

                // initialize the custom_index property values
                graph_traits<Graph>::vertex_iterator vi, vend;
                graph_traits<Graph>::vertices_size_type cnt = 0;

                for(boost::tie(vi,vend) = vertices(read_g); vi != vend; ++vi) {
                    put(index, *vi, cnt++);
                }

                /*
                
                END OF READ
                
                */

                //std::vector<double> tol = {0.1, 0.2, 0.3};
                double tol = 0.2; // for testing purposes
                bool debug = false;
                // calculate QSP BWss

                auto start = high_resolution_clock::now();
                    
                compute_tolSP_paths(read_g, tol, graph_number, debug);
                
                auto stop = high_resolution_clock::now();
                completed++;
                std::chrono::duration<float> duration = (stop - start);
                std::cout <<"Time taken for this graph: "
                    << duration.count() << " seconds" << std::endl;

                std::cout << "Completed " << completed <<  " graphs" << std::endl;
                char out_char[50];
                std::sprintf(out_char, "Completed %d graphs\n", completed);
                log(out_char);
                //break; // Remove this line to process all files in the directory
            }
        }
    }
    return 0;
}

