#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <fstream>
#include <map>
#include <chrono>
#include <numeric>
#include <future>

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



// #define FMT_DEPRECATED_OSTREAM
// #include <fmt/ranges.h>
// #include <fmt/ostream.h>



using namespace boost;
using namespace std::chrono;

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

using Filtered = filtered_graph<Graph, boost::function<bool(Edge)>, boost::function<bool(Vertex)> >;



// concatenates root and spur 
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

// FIND SPUR PATH FUNCITON
std::tuple<Path, double> FindShortestPathYen(Graph &G, int spurNode, int target,
                                                            const std::vector<EdgePair> &devLinks,
                                                            const Path &rootPath){
// Graph GOrig = G;

// have a filtered "copy" f that just removes a set of vertices:
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

// std::unordered_map<int, Graph::vertex_descriptor> indexToVertexDescriptorMap;
// graph_traits<Graph>::vertex_iterator vi, vend;
// for (tie(vi, vend) = vertices(G); vi != vend; ++vi) {
//     int index = get(&VertexProperties::custom_index, G[*vi]); // Get the index property from VertexProperties
//     indexToVertexDescriptorMap[index] = *vi;
// }

// Vertex spurVertex = indexToVertexDescriptorMap[spurNode];

int spurDeg = degree(spurNode, f);
        if (spurDeg == 0){
            return std::make_tuple(spurPath, std::numeric_limits<double>::max());
        };

// Vertex targVertex = indexToVertexDescriptorMap[target];

int targDeg = degree(target, f);
if (targDeg == 0){
    return std::make_tuple(spurPath, std::numeric_limits<double>::max());
};

// Dijkstra's algorithm to find the spur path

std::size_t numVertices = std::distance(boost::vertices(G).first, boost::vertices(G).second);

Path predecessors(numVertices);
std::vector<double> distances(numVertices);
        
dijkstra_shortest_paths(f, spurNode,
                        predecessor_map(make_iterator_property_map(predecessors.begin(), get(vertex_index, f)))
                        .distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, f))));

// Reconstruct the spur path

int currentNode = target;

while (currentNode != spurNode) {

    spurPath.push_back(currentNode);
    if (predecessors[currentNode] == currentNode) {
        spurPath.clear();
        return std::make_tuple(spurPath, std::numeric_limits<double>::max());
    }
    currentNode = predecessors[currentNode];
}

spurPath.push_back(spurNode);
std::reverse(spurPath.begin(), spurPath.end());
removed_set.clear();
removed_edge_set.clear();
return std::make_tuple(spurPath, distances[target]);

}

// Finds deviation links given previous paths and root path
std::vector<EdgePair> findDeviationLinks(Graph G, std::vector<Path> paths, Path rootPath, int sizeOfRoot){

    std::vector<EdgePair> devLinks;
    for (Path ipath : paths) {
        if (ipath.size() > sizeOfRoot + 1){
            Path prev_path(ipath.begin(), ipath.begin() + sizeOfRoot + 1);

            if (rootPath == prev_path){
                // printf("before findspur\n");
                //finding deviation links
                Vertex u = ipath[sizeOfRoot];
                Vertex v = ipath[sizeOfRoot + 1];

                //check if dev. link is in graph
                Edge e;
                bool exists = false;
                std::pair<Edge, bool> edgeInfo = edge(u, v, G);
                e = edgeInfo.first;
                exists = edgeInfo.second;

                if (exists) {
                    devLinks.emplace_back(u, v);
                } else {
                    printf("deviation link (%lu, %lu) not in graph\n", v, u);
                }
                // printf("after findspur\n");
            }
        }
    }
    return devLinks;
}

/*

******************** YEN'S MAIN FUNCTION KSP ************************

*/ 

std::tuple<std::vector<Path>, std::vector<double>> YenTolShortestPaths(Graph &G, Vertex source, Vertex target, const double tol) {
std::priority_queue<std::pair<double, Path >,
std::vector<std::pair<double, Path > >,
std::greater<std::pair<double, Path > > > B; //heap of path weight, path tuples

std::vector<double> pathLens;
std::vector<Path> paths;

std::vector<double> lengths = std::vector<double>(num_vertices(G), std::numeric_limits<double>::infinity());//vector of SP from dijkstra

paths.clear();
lengths[source] = 0;//set length to source to 0

Path path(num_vertices(G), -1);  // Initialize path with -1 for all vertices
std::vector<double> distances(num_vertices(G));  // Store distances for Dijkstra's algorithm

dijkstra_shortest_paths(G, source,
predecessor_map(make_iterator_property_map(path.begin(), get(vertex_index, G)))
.distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, G))));

// Extract the path to the target node
Path pathToTarget;
for (Vertex v = target; v != source; v = path[v]) {
    pathToTarget.push_back(v);
}
pathToTarget.push_back(source);
double dist_SP = distances[target];//distance of the shortest path

// Reverse the path to get it from source to target
std::reverse(pathToTarget.begin(), pathToTarget.end());

paths.push_back(pathToTarget);
pathLens.push_back(dist_SP);
lengths[target] = dist_SP;

struct VectorHash {
std::size_t operator()(const Path& vec) const {
std::size_t seed = vec.size();
for(const auto& i : vec) {
    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
return seed;
}
};

std::unordered_set<Path, VectorHash> elementsInQueue;//to check for repeated paths in the heap


double SP = lengths[target];
double l = SP;

int count = 0;
while (l <= SP * (1 + tol)) {
    // printf("%lf\n", l);
    count++;
    
    for (int j = 0; j < paths.back().size() - 1; ++j) {
        Vertex spurNode = paths.back()[j];
        Path rootPath(paths.back().begin(), paths.back().begin() + j + 1);//takes the root path from the last path using iterators

        //iterates over all paths

        std::vector<EdgePair> devLinks = findDeviationLinks(G, paths, rootPath, j);

        int spurDeg = degree(spurNode, G);
        if (spurDeg == devLinks.size()){
            continue;
        };
        

        std::tuple<Path, double> spurPath_t = FindShortestPathYen(G, spurNode, target, 
                                                                                devLinks, rootPath);                          

        Path spurPath = std::get<0>(spurPath_t);
        double spurPath_len = std::get<1>(spurPath_t);

        if (!spurPath.empty() && spurPath.back() == target) 
            {

            //concatenating and retrieving total path and total path length
            std::tuple totalPath_t = concatenatePaths(G, rootPath, spurPath); 

            Path totalPath = std::get<0>(totalPath_t);
            double totalPathLength = std::get<1>(totalPath_t);


            if (elementsInQueue.find(totalPath) == elementsInQueue.end()) {
                elementsInQueue.insert(totalPath);
                B.push(std::make_pair(totalPathLength, totalPath));
            }

        }
    }
    
    if (!B.empty()) {
        auto lp = B.top();
        B.pop();
        l = lp.first;
        auto p = lp.second;

        if (l <= SP * (1 + tol)) {
            // printf("closeness to Tolerance = %lf\n", (l-SP)/(SP * tol));
            lengths[target] = l;
            paths.push_back(p);
            pathLens.push_back(l);
        }
    } else {
        break;
    }
}
return std::make_tuple(paths, pathLens);
}


std::vector<double> compute_pair(Graph G, Vertex n1, Vertex n2, double tol, Path ini_pred, std::vector<double> ini_dists) {
    int nnodes = num_vertices(G);
    std::vector<double> betw_pair_1(nnodes, 0.0);

    // Populate node_list with node indices
    Path node_list(nnodes); 
    Graph::vertex_iterator vi, vend;
    for (std::tie(vi, vend) = vertices(G); vi != vend; ++vi){
        node_list[*vi] = *vi;
    }
    
    auto pair_t = YenTolShortestPaths(G, n1, n2, tol);

    auto paths = std::get<0>(pair_t);
    int num_paths = paths.size();

    for (const auto& p : paths) {
        for (int i = 1; i < p.size() - 1; ++i) {
            auto it = std::find(node_list.begin(), node_list.end(), p[i]);
            if (it != node_list.end()) {
                int index = std::distance(node_list.begin(), it);
                betw_pair_1[index] += 1.0;
            }
        }
    }

    for (int i = 0; i < node_list.size(); ++i) {
        betw_pair_1[i] /= double(num_paths);
    }
    printf("number of paths: %zu, pair: %lu, %lu\n", paths.size(), n1, n2);
    return betw_pair_1;
}
template<typename T>
bool isReady(const std::future<T>& f) {
        if (f.valid()) { // otherwise you might get an exception (std::future_error: No associated state)
            return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        } else {
            return false;
        }
    }

std::vector<double> compute_tolSP_betw(Graph &G, double tol) {
    int nnodes = num_vertices(G);
    std::vector<double> betw(nnodes, 0);
    std::vector<std::future<std::vector<double>>> futures;
    std::vector<std::unique_ptr<Graph>> graphs; // Store unique pointers for each graph
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

    it = 0;
    #pragma omp parallel for private(ini_path_t) shared(betw, it)
    for (int k = 0; nnodes*(nnodes-1)/2 > k; k++) {
        auto [i, j] = loop_map[k];
        ini_path_t = dijkstra_dict.at(j);
        std::vector<double>  result = compute_pair(G, i, j, tol, std::get<0>(ini_path_t), std::get<1>(ini_path_t));

        // This critical section ensures that only one thread at a time can update betw.
        #pragma omp critical
        {
            it++;
            printf("%d, %d\n", i, j);
            std::transform(betw.begin(), betw.end(), result.begin(), betw.begin(), std::plus<double>());
            printf("progress: %lf\n", 2*it/double(nnodes*(nnodes-1)));
            }
        
    }

    #pragma omp barrier
    printf("tasks completed\n");
    return betw;
}


int main() {

    /*
    
    READING CSV FILE: EDGE[0], EDGE[1], WEIGHT

    */

    Graph delaunay(0);

    std::ifstream file("time_testing_delaunays/delaunay300.csv");
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        boost::tokenizer<boost::escaped_list_separator<char>> tokens(line);
        auto tokenIterator = tokens.begin();
        int vertex1 = std::stoi(*tokenIterator++);
        int vertex2 = std::stoi(*tokenIterator++);
        double weight = std::stod(*tokenIterator);
        // Add edge to the graph with the given weight
        Edge e = add_edge(vertex1, vertex2, delaunay).first;
        put(edge_weight, delaunay, e, weight);
        // boost::add_edge(vertex1, vertex2, EdgeProperties{weight}, delaunay);
    }

    /*
    
    END OF READ
    
    */

    // Find K shortest paths 

    // std::vector<double> times(100);

    // // for (int i = 0; 100>i; i++){

        auto start = high_resolution_clock::now();

        Path ini_p;
        std::vector<double> dist_ini;
        std::pair<std::vector<Path>, std::vector<double>> pair_t = YenTolShortestPaths(delaunay, 0, 290, 0.1);
        

        std::vector<Path> paths = std::get<0>(pair_t);
        std::vector<double> pathLens = std::get<1>(pair_t);
        auto stop = high_resolution_clock::now();

        std::chrono::duration<float> duration = (stop - start);

    //     // times[i] = duration.count();
        printf("number of paths: %zu\n", paths.size());

        std::cout << "Time taken by function: "
            << duration.count() << " seconds" << std::endl;

    // // }

    // double mean_t = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    // double sq_sum = std::inner_product(times.begin(), times.end(), times.begin(), 0.0);
    // double stdev = std::sqrt(sq_sum / times.size() - mean_t * mean_t);
    // printf("mean time: %lf ± %lf", mean_t, stdev);

    // Output the results

    // for (int j = 0; j < yen.paths.size();++j) {
    //     const auto &path = yen.paths[j];
    //     std::cout << "Path Length: " << yen.pathLens[j] << ", Path: ";
    //     for (int i = 0; i < path.size(); ++i) {
    //         std::cout << path[i];
    //         if (i < path.size() - 1) {
    //             std::cout << " -> ";
    //         }
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "Time taken by function: "
    //     << duration.count() << " seconds" << std::endl;
    

    // auto start = high_resolution_clock::now();
    // double tol = 0.1;
    // auto betweenness = compute_tolSP_betw(delaunay, tol);
    
    // auto stop = high_resolution_clock::now();

    // std::chrono::duration<float> duration = (stop - start);
    // // times[i] = duration.count();
    
    // // file pointer 
    // std::fstream fout; 
  
    // // opens an existing csv file or creates a new file. 
    // fout.open("tol_betw_VF.csv", std::ios::out); 

    // std::cout << "Tol Betweenness: \n";
    // for (int i = 0; i < betweenness.size(); ++i){
    //     std::cout << "node " << i << ": " << betweenness[i] << " ";
    //     fout << i << "," << betweenness[i] << "\n";
    // }
    // std::cout << "\n"; 

    // std::cout << "Time taken by function: "
    //    << duration.count() << " seconds" << std::endl;
    return 0;
}
    





