#ifndef QSP_FUNCTIONS_H
#define QSP_FUNCTIONS_H

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <type_traits>
#include <vector>
#include <functional>
#include <omp.h>

#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/process.hpp>

using namespace boost;

struct VertexProperties {
    int custom_index;
};

typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    VertexProperties,
    boost::property<boost::edge_weight_t, double>
> Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef std::pair<Vertex, Vertex> EdgePair;
typedef std::vector<Vertex> Path;

struct VertexProps {
    Path path;
    int id;
};

typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS,
    VertexProps
> Tree;

typedef graph_traits<Tree>::vertex_descriptor VertexT;
typedef graph_traits<Tree>::edge_descriptor EdgeT;
typedef std::unordered_map<EdgePair, Tree, boost::hash<EdgePair>> TreeDict;

using Filtered = filtered_graph<Graph, boost::function<bool(Edge)>, boost::function<bool(Vertex)> >;

class Heuristic :
      public boost::astar_heuristic<Filtered, double>
{
public:
    explicit Heuristic(const std::vector<double> distances) : dist_vec(distances) {}
    double operator()(Vertex v) { return dist_vec[v]; }
private:
    const std::vector<double> dist_vec;
};

struct found_goal {};

class astar_goal_visitor : public boost::default_astar_visitor
{
public:
    explicit astar_goal_visitor(Vertex goal) : m_goal(goal) {}
    template <class FilteredGraph>
    void examine_vertex(Vertex u, FilteredGraph&)
    {
        if (u == m_goal){
            throw found_goal();
        }
    }
private:
    Vertex m_goal;
};

// struct to find all predecessors in unitary weights using dijkstra
struct state {
            Graph&                       g;
            std::vector<double>      dist;
            std::vector<std::set<Vertex>> pred;
            state(Graph& g) : g(g), dist(num_vertices(g)), pred(num_vertices(g)) {}
        };

// Dijkstra visitor to get all shortest paths
struct vis_t : boost::default_dijkstra_visitor {
    state& s;
    vis_t(state& s) : s(s) {}

    void edge_relaxed(Edge e, Graph const& g) const {
        s.pred[target(e, g)] = {source(e, g)};
    }
    void edge_not_relaxed(Edge e, Graph const& g) const {
        // e: u -> v
        auto u = source(e, g);
        auto v = target(e, g);

        auto old  = s.dist[v];
        auto new_ = s.dist[u] + get(boost::edge_weight, g, e);

        if (std::nextafter(new_, old) == old)
            s.pred[v].insert(u);
    }
}; 


/*--------------- FUNCTIONS FOR SPR ALGORITHM ---------------*/

/**
 * @brief Reconstruct every shortest path from a predecessor-set table.
 *
 * Performs a depth-first enumeration (using an explicit stack) to walk the
 * predecessor sets recorded during multi-predecessor Dijkstra, returning all
 * distinct shortest paths between the supplied endpoints.
 *
 * @param predecessors  Vector of predecessor sets per vertex.
 * @param initialNode   Source vertex.
 * @param destinationNode Destination vertex.
 * @return Collection of source-to-destination paths expressed as vertex IDs.
 */
std::vector<std::vector<Vertex>> getPaths(const std::vector<std::set<Vertex>>& predecessors, 
                                        Vertex initialNode, Vertex destinationNode) {
    std::vector<std::vector<Vertex>> paths;
    std::stack<std::vector<Vertex>> stack;

    stack.push({destinationNode});

    while (!stack.empty()) {
        std::vector<Vertex> currentPath = stack.top();
        stack.pop();

        int currentNode = currentPath.back();

        if (currentNode == initialNode) {
            // Found a path
            std::reverse(currentPath.begin(), currentPath.end());
            paths.push_back(currentPath);
        } else {
            // Backtrack using predecessors
            for (int predecessor : predecessors[currentNode]) {
                std::vector<Vertex> newPath = currentPath;
                newPath.push_back(predecessor);
                stack.push(newPath);
            }
        }
    }

    return paths;
}

/**
 * @brief Concatenate a root path and spur path, returning the combined length.
 */
inline std::tuple<Path,double> concatenatePaths(Graph G, const Path &path1, const Path &path2) {

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

/**
 * @brief Run the A* spur search used by Yen's algorithm on a filtered graph.
 */
inline std::tuple<Path, double> FindShortestPathYen(Graph &G, int spurNode, int target,
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
            
    // dijkstra_shortest_paths(f, spurNode,
    //                         predecessor_map(make_iterator_property_map(predecessors.begin(), get(&VertexProperties::custom_index, f)))
    //                         .distance_map(make_iterator_property_map(distances.begin(), get(&VertexProperties::custom_index, f))));

    // // Reconstruct the spur path

    // int currentNode = target;

    // while (currentNode != spurNode) {

    //     spurPath.push_back(currentNode);
    //     if (predecessors[currentNode] == currentNode) {
    //         spurPath.clear();
    //         return std::make_tuple(spurPath, std::numeric_limits<double>::max());
    //     }
    //     currentNode = predecessors[currentNode];
    // }
    // spurPath.push_back(spurNode);
    // std::reverse(spurPath.begin(), spurPath.end());

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

/**
 * @brief Retrieve a spur candidate stored in the initial shortest-path tree.
 */
inline Path GetCandidatePathVF(const Graph& G, std::vector<EdgePair> dev_links, 
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

/**
 * @brief Retrieve a valid spur path from the candidate tree, if one exists.
 */
inline std::tuple<Path, Path, VertexT> RetrieveCandidatePath(Tree& candTree, const std::vector<Vertex>& N, Vertex devNode) {
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


/**
 * @brief Expand the candidate tree to grow a spur path when reuse fails.
 */
inline Path CalculateCandidatePaths(Graph& G, Tree& candTree, const Path& N,
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
            // printf("curr node %d, num nodes %d\n", current_node, num_vertices(candTree));

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
    // dijkstra_shortest_paths(f, s,
    //                         predecessor_map(make_iterator_property_map(p.begin(), get(&VertexProperties::custom_index, f)))
    //                         .distance_map(make_iterator_property_map(d.begin(), get(&VertexProperties::custom_index, f))));
    // // If deviation node is in the candidate path and no loop is formed, update the candidate path
    // candPath.clear();
    // Vertex prev_v = graph_traits<Graph>::null_vertex();
    // for (Vertex v = target; v != s; v = p[v]) {
    //     candPath.push_back(static_cast<Vertex>(v));
    //     if (prev_v == v){
    //         Path empt;
    //         return empt;
    //     }
    //     prev_v = v;
    // }
    // // candPath.push_back(static_cast<int>(s));
    // reverse(candPath.begin(), candPath.end());   

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

/**
 * @brief Find or recompute the spur path for the current deviation node.
 */
inline std::pair<Path, int> FindSpurPath_SPR(Graph& G, const std::vector<EdgePair>& devLinks, const Path& N, Vertex devNode,
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
        } catch (boost_swap_impl::out_of_range&) {
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

/**
 * @brief Gather deviation edges for all paths sharing the supplied root.
 */
inline std::vector<EdgePair> findDeviationLinks(Graph G, std::vector<Path> paths, Path rootPath, int sizeOfRoot){

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
                } else {
                    // printf("deviation link (%lu, %lu) not in graph\n", v, u);
                }
                // printf("after findspur\n");
            }
        }
    }
    return devLinks;
}

/*

************************ SPR MAIN FUNCTION **************************

*/

/**
 * @brief Compute tolerance-respecting paths using predecessor sets (unit/weighted graphs).
 *
 * This routine enumerates every shortest path discovered via the predecessor-set
 * Dijkstra visitor, then continues Yen's enumeration within the tolerance band.
 *
 * @param G                 Working graph.
 * @param source            Source vertex.
 * @param target            Target vertex.
 * @param ini_predecessors  Predecessor sets (populated if empty).
 * @param distances         Distance vector (populated if empty).
 * @param K                 Number of K-shortest paths.
 * @return Tuple holding tolerance-respecting paths (with original IDs) and their lengths.
 */
inline std::tuple<std::vector<Path>, std::vector<double>> get_KSPs(Graph G, Vertex source, Vertex target, std::vector<std::set<Vertex>>& ini_predecessors, 
                                                                        std::vector<double>& distances, const int K_max) {
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
        state s1(G);
        vis_t vis{s1};

        dijkstra_shortest_paths(G, target,
                                weight_map(get(edge_weight, G))
                                .distance_map(make_iterator_property_map(s1.dist.begin(), get(vertex_index, G)))
                                .visitor(vis));

        ini_predecessors.insert(ini_predecessors.end(), s1.pred.begin(), s1.pred.end());
        distances.insert(distances.end(), s1.dist.begin(), s1.dist.end());
    }

    // Extract the path to the target node
    
    std::vector<Path> pathsToTarget = getPaths(ini_predecessors, target, source);
    double dist_SP = distances[source]; //distance of the shortest path
    
    // -----

for (Path ini_p: pathsToTarget){
    std::reverse(ini_p.begin(), ini_p.end());
    real_paths.emplace_back(ini_p);
    pathLens.emplace_back(dist_SP);    
}

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

    // // Iterate through the vector using a range-based for loop
    // for (int i = 0; i < dists_to_s.size(); ++i) {
    //     if ((distances[i] + dists_to_s[i]) > dist_SP*(1+tol)) {
    //         outranged_n.insert(i);
    //     }
    // }

    Graph reduced_G;
    copy_graph(G_red, reduced_G);

    printf("reduced total number of nodes from %zu to %zu\n", num_vertices(G), num_vertices(reduced_G));

    // Access the vertex index property map
    auto vertexIndexMap = get(boost::vertex_index, reduced_G);

    // map of old ID to new ID
    std::unordered_map<std::size_t, Vertex> inverseMap;

    // Iterate over vertices and print node_id property
    boost::graph_traits<Graph>::vertex_iterator vertexIt, vertexEnd;
    for (boost::tie(vertexIt, vertexEnd) = boost::vertices(reduced_G); vertexIt != vertexEnd; ++vertexIt) {
        std::size_t vertexID = boost::get(vertexIndexMap, *vertexIt);
        std::size_t nodeID = reduced_G[*vertexIt].custom_index;
        inverseMap[nodeID] = vertexID;
        // std::cout << "Vertex " << vertexID << ": Node ID = " << nodeID << std::endl;
    }

    // New source and target IDs

    Vertex new_s = inverseMap[source];
    Vertex new_t = inverseMap[target];

    //Do the initial dijkstra with the new node ids

    Path ini_preds;
    std::vector<double> dists;
    ini_preds.resize(num_vertices(reduced_G), -1);  // Initialize path with -1 for all vertices
    dists.resize(num_vertices(reduced_G), 0); // Store distances for Dijkstra's algorithm

    dijkstra_shortest_paths(reduced_G, new_t,
                            predecessor_map(make_iterator_property_map(ini_preds.begin(), get(vertex_index, reduced_G)))
                            .distance_map(make_iterator_property_map(dists.begin(), get(vertex_index, reduced_G))));

    

    // Extract the path to the new target node id
    Path new_SP;

    for (Vertex v = new_s; v != new_t; v = ini_preds[v]) {
        
        new_SP.emplace_back(v);

    }
    new_SP.emplace_back(new_t);
    // -----

    paths.emplace_back(new_SP);

    double l = dist_SP;
    int tot_spur = 0;
    int n_reused = 0;

    for (int k = 1; k < K_max; k++){
    // while (l <= dist_SP * (1 + tol)) {

        for (int j = 0; j < paths.back().size() - 1; ++j) {
            Vertex spurNode = paths.back()[j];
            Path rootPath(paths.back().begin(), paths.back().begin() + j + 1);//takes the root path from the last path using iterators

            //iterates over all paths

            std::vector<EdgePair> devLinks = findDeviationLinks(reduced_G, paths, rootPath, j);

            // if the degree of the spur node is the same as the amount of deviation 
            // links continue with the next spur node
            int spurDeg = degree(spurNode, reduced_G);
            
            if (spurNode == source){
                
                if (spurDeg == devLinks.size()){
                    continue;
                };
                
            } else {

                if (devLinks.size() >= (spurDeg-1)){
                    continue;
                };
            }

            
            Path N(paths.back().begin(), paths.back().begin() + j);
            tot_spur++;

            auto spurPath_tp = FindSpurPath_SPR(reduced_G, devLinks, N, spurNode, dists, 
                                            ini_preds, treeMap, rootPath, new_s, new_t);  
            

            Path spurPath = std::get<0>(spurPath_tp);
            n_reused += std::get<1>(spurPath_tp);
            
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
                // if (totalPathLength <= dist_SP * (1 + tol)){

                //finding if total path is in the heap and adding it accordingly.
                auto it = std::find_if(B.begin(), B.end(), [&totalPath](const auto& element) {
                return std::get<1>(element) == totalPath;
                });
                spurPath.clear();

                if (it == B.end() || B.size() == 0) {
                    B.emplace_back(std::make_pair(totalPathLength, totalPath));
                    std::push_heap(B.begin(), B.end(), CompareTuples());

                // }
                }
            }
        }


        if (!B.empty()) {
            auto lp = B.front();
            std::pop_heap(B.begin(), B.end(), CompareTuples());
            B.pop_back();
            l = std::get<0>(lp);
            auto p = std::get<1>(lp);

            // for (int i = 0; i < p.size(); ++i) {
            //     std::cout << p[i] << " ";
            //     if (i < p.size() - 1) {
            //         std::cout << "-> ";
            //     }
        
            // }

            // for (int i = 0; i < p.size()-1; ++i) {
            //     Edge edge = boost::edge(p[i], p[i+1], G).first;

            //     std::cout << "(" << p[i] << ", " << p[i+1] << ")" << ": " << boost::get(boost::edge_weight, G, edge) << "    ";

            // }
            // printf(" l = %lf\n", l);

            // if (l <= dist_SP * (1 + tol)) {
            //     printf("closeness to Tolerance = %lf\n", (l-dist_SP)/(dist_SP * tol));
                // printf("k = %d\n", k);
            paths.emplace_back(p);

            Path real_p;
            for (int i = 0; i < p.size(); ++i) {
                real_p.emplace_back(reduced_G[p[i]].custom_index);
            }
            real_paths.emplace_back(real_p);
            pathLens.emplace_back(l);

            // }
        } else {
            break;
        }
    }
    //printf("reused ratio %lf\n", n_reused/double(tot_spur));
    return std::make_tuple(real_paths, pathLens);
}


/**
 * @brief Compatibility overload that mirrors the original predecessor-vector API.
 */
inline std::tuple<std::vector<Path>, std::vector<double>> TolShortestPaths_SPR(Graph G, Vertex source, Vertex target, Path& ini_predecessors, 
                                                                        std::vector<double>& distances, const double tol) {
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

    // std::cout << std::endl;
    Graph reduced_G;
    copy_graph(G_red, reduced_G);

    // printf("reduced total number of nodes from %zu to %zu\n", num_vertices(G), num_vertices(reduced_G));

    // Access the vertex index property map
    auto vertexIndexMap = get(boost::vertex_index, reduced_G);

    // map of old ID to new ID
    std::unordered_map<std::size_t, Vertex> inverseMap;

    // Iterate over vertices and print node_id property
    boost::graph_traits<Graph>::vertex_iterator vertexIt, vertexEnd;
    for (boost::tie(vertexIt, vertexEnd) = boost::vertices(reduced_G); vertexIt != vertexEnd; ++vertexIt) {
        std::size_t vertexID = boost::get(vertexIndexMap, *vertexIt);
        std::size_t nodeID = reduced_G[*vertexIt].custom_index;
        inverseMap[nodeID] = vertexID;
        // std::cout << "Vertex " << vertexID << ": Node ID = " << nodeID << std::endl;
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

    

    // Extract the path to the new target node id
    Path new_SP;

    for (Vertex v = new_s; v != new_t; v = ini_preds[v]) {
       //printf("%lu, %lu\n", v, new_t);
        
        new_SP.emplace_back(v);

    }
    new_SP.emplace_back(new_t);
    // -----

    paths.emplace_back(new_SP);

    double l = dist_SP;
    // int tot_spur = 0;
    // int n_reused = 0;

    // for (int k = 1; k < 2000; k++){
    while (l <= dist_SP * (1 + tol)) {

        for (int j = 0; j < paths.back().size() - 1; ++j) {
            Vertex spurNode = paths.back()[j];

            Path rootPath(paths.back().begin(), paths.back().begin() + j + 1);//takes the root path from the last path using iterators

            //iterates over all paths

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

            
            Path N(paths.back().begin(), paths.back().begin() + j);
            // tot_spur++;

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

            // for (int i = 0; i < p.size(); ++i) {
            //     std::cout << p[i] << " ";
            //     if (i < p.size() - 1) {
            //         std::cout << "-> ";
            //     }
        
            // }

            // for (int i = 0; i < p.size()-1; ++i) {
            //     Edge edge = boost::edge(p[i], p[i+1], G).first;

            //     std::cout << "(" << p[i] << ", " << p[i+1] << ")" << ": " << boost::get(boost::edge_weight, G, edge) << "    ";

            // }
            // printf(" l = %lf\n", l);

            if (l <= dist_SP * (1 + tol)) {
                // printf("closeness to Tolerance = %lf\n", (l-dist_SP)/(dist_SP * tol));
                // printf("k = %d\n", k);
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
    }
    // printf("reused ratio %lf\n", n_reused/double(tot_spur));
    return std::make_tuple(real_paths, pathLens);
}




/*------------------------------------TOLERANCE BETWEENNESS------------------------------------*/

/**
 * @brief Compute tolerance betweenness contributions for a single source-target pair.
 *
 * @tparam doub_vec  Either a scalar tolerance or a vector of tolerances.
 */
template<typename doub_vec>
inline std::vector<std::vector<double>> compute_pair(Graph G, Vertex n1, Vertex n2, doub_vec tols, Path ini_pred, std::vector<double> ini_dists) {

    // choose the maximum tol among the vector
    double tol = 0;
    bool isvec = true;
    int it = 0;
    int size_tols = 1;
    std::vector<double> tols_tot(size_tols, 0.0);

    if constexpr (std::is_same_v<doub_vec, std::vector<double>>) {

        isvec = true;
        tol = *std::max_element(tols.begin(), tols.end());
        // if 0% tolerance is not in tols we add it
        if (std::find(tols.begin(), tols.end(), 0) == tols.end()) {
            size_tols = tols.size()+1;     
            tols_tot.resize(size_tols, 0.0);
            std::copy(tols.begin(), tols.end(), tols_tot.begin() + 1);
        } else {
            size_tols = tols.size();
            tols_tot.resize(size_tols, 0.0);
            std::copy(tols.begin(), tols.end(), tols_tot.begin());
            }
    } else if constexpr (std::is_same_v<doub_vec, double>) {

        isvec = false;
        tol = tols;

    } else {

        printf("tolerance type is not valid\n");
        exit(0);
    }

    int nnodes = num_vertices(G);
    // the betweenness contribution vec will always contain the 0% tolerance and then all the tols specified in the vector
    std::vector<std::vector<double>> betw_pair(size_tols, std::vector<double>(nnodes, 0));
    
    // compute the tolerance KSPs
    std::vector<double> ini_dists2 = ini_dists;
    Path ini_pred2 = ini_pred;
    auto pair_t = TolShortestPaths_SPR(G, n1, n2, ini_pred, ini_dists, tol);

    std::vector<double> dists = std::get<1>(pair_t);

    auto paths = std::get<0>(pair_t);
    int num_paths = paths.size();

    double d_sp = dists[0];

    // add the BWss contributions to each node

    int low_tol_ind = 0;
    int new_ind = 0;
    std::vector<double> n_tol_ksps(size_tols, 0); // number of K Shortest Paths in each tolerance

    // the number of paths of the largest tolerance is the total number of recorded paths
    n_tol_ksps[size_tols-1] = num_paths;

    for (const auto& p : paths) {

        double d = dists[it];
        // for each path we find the lowest tolerance that does not contain the path length
        // First we update the number of shortest paths if the distance is still the SP
        // if (d == d_sp){
        //     n_tol_ksps[0] = it+1;
        // }

        if (isvec){
            // As we iterate through the KSPs, we check if the path length falls below the input tolerances
            for (int i = low_tol_ind; i<size_tols; i++){
                // We find the lowest tolerance that contains the path length
                if (d <= d_sp*(1+tols_tot[i])){
                    // when we exceed one of the tolerances we update the number of KSP that 
                    // fall in that tolerance to later normalize each tol BW
                    
                    // if (i == 0){
                    //     printf("d-1 %lf, d %lf, prc = %lf, n_p %d, tot_p %d\n", dists[it-1], d, d_sp*1.1, it, num_paths);
                    // }

                    n_tol_ksps[i] = it+1;

                    //if the next path jumps one or more of the tolerances (eg. from 0 to 30%), 
                    //update also the number of paths for the previous tolerances
                    // CHECK!!!!
                    // for (int j = 1; j<i+1; j++){
                    //     if (n_tol_ksps[j] == 0){
                    //         n_tol_ksps[j] = it;
                    //     }
                    // }
                } else {
                    // if the path length is larger than the current tolerance, we update the index of the lowest tolerance
                    new_ind = i+1;
                }

            }
            // update the index of the lowest tolerance 
            low_tol_ind = new_ind;
        }

        for (int i = 1; i < p.size() - 1; ++i) {

            // // fill up the SP BW only if the path has the shortest length
            // if (d == d_sp){
            //     betw_pair[0][p[i]] += 1.0;
            // }

            // update the contributions of the BW for different tolerances
            for (int k = low_tol_ind; k < size_tols; k++){
                betw_pair[k][p[i]] += 1.0;
                
            }
        }
        it++;
    }
    // for (int j = 1; j<size_tols+1; j++){
    //     if (n_tol_ksps[j] == 0){
    //         n_tol_ksps[j] = num_paths;
    //     }
    // }
    for (int k = 0; k < size_tols; k++){
        for (int i = 0; i < nnodes; ++i) {
            betw_pair[k][i] /= n_tol_ksps[k];
        }
}
    return betw_pair;
}

// tol SP Betw function
template<typename doub_vec>

/**
 * @brief Aggregate tolerance betweenness values for every vertex in the graph.
 */
inline std::vector<std::vector<double>> compute_tolSP_betw(Graph &G, doub_vec tols) {

    int nnodes = num_vertices(G);
    int size_tols = 1;
    if constexpr (std::is_same_v<doub_vec, std::vector<double>>) {
        size_tols = tols.size();
        // sort the tolerances from small to large for later
        std::sort(tols.begin(), tols.end());
        
    } else if constexpr (std::is_same_v<doub_vec, double>) {

        size_tols = 1;

    } else {

        printf("tolerance type is not valid\n");
        exit(0);
    }

    // array of contributions to a certain node or group of nodes (ring)

    // std::vector<std::vector<double>> pair_contr_arr(nnodes, std::vector<double>(nnodes, 0));

    std::vector<std::vector<double>> betw(tols.size()+1, std::vector<double>(nnodes, 0));
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
    #pragma omp parallel for private(ini_path_t) shared(betw, it, tols)
    for (int k = 0; nnodes*(nnodes-1)/2 > k; k++) {
        auto [i, j] = loop_map[k];
        ini_path_t = dijkstra_dict.at(j);
        std::vector<std::vector<double>> result = compute_pair(G, i, j, tols, std::get<0>(ini_path_t), std::get<1>(ini_path_t));

        // This critical section ensures that only one thread at a time can update betw.
        #pragma omp critical
        {
            it++;
            // printf("%d, %d\n", i, j);
            // double max_contibution = 0;

            // for (int n = 0; n<8 ; n++){

            //     if (result[n] > max_contibution) {

            //         max_contibution += result[n];
            //     }
            // }
            // pair_contr_arr[i][j] = max_contibution;
            // pair_contr_arr[j][i] = max_contibution;

            for (int tol_ind = 0; tol_ind < size_tols+1; tol_ind++){
                std::transform(betw[tol_ind].begin(), betw[tol_ind].end(), result[tol_ind].begin(), betw[tol_ind].begin(), std::plus<double>());
            }
            // printf("progress: %lf\n", 2*it/double(nnodes*(nnodes-1)));
            }
        
    }

    #pragma omp barrier
    printf("tasks completed\n");
    return betw;
}

/**
 * @brief Extract the first numeric substring from a file name.
 */
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

/**
 * @brief Write betweenness results for each tolerance into CSV files.
 */
inline void writeBTWtoCSV(const std::vector<std::vector<double>>& data, const std::string& rootPath, std::vector<double> tols) {
    // Ensure the directory exists
    int size_tols = 1;
    std::vector<double> tols_tot(size_tols, 0.0);
    // if 0% tolerance is not in tols we add it
    if (std::find(tols.begin(), tols.end(), 0) == tols.end()) {
        size_tols = tols.size()+1;
        tols_tot.resize(size_tols, 0.0);
        std::copy(tols.begin(), tols.end(), tols_tot.begin()+1);
    } else {
        size_tols = tols.size();
        std::vector<double> tols_tot = tols;
        }

    for (size_t i = 0; i < data.size(); ++i) {
        printf("Writing data for tolerance %.1f\n", tols_tot[i]*100);
        // Construct full path: rootPath/vector_i.csv
        std::string filename = rootPath + "/qspbw_tol_" + std::to_string(static_cast<int>(tols_tot[i]*100)) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            continue;
        }

        // Write header
        file << "index,value\n";

        // Write data
        const auto& vec = data[i];
        for (size_t j = 0; j < vec.size(); ++j) {
            file << j << "," << vec[j] << "\n";
        }

        file.close();
    }
}

#endif // QSP_FUNCTIONS_H
