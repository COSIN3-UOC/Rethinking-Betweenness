#include "qsp_bw.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/property_map/property_map.hpp>
#include <random>
#include <utility>
#include <boost/tokenizer.hpp>
#include <omp.h>


using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return (static_cast<size_t>(p.first) << 32) ^ static_cast<size_t>(p.second);
    }
};

// ---------- all-pairs shortest paths ----------
std::vector<std::vector<double>> all_pairs_shortest_paths(const Graph& g) {
    const size_t n = boost::num_vertices(g);
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> D(n, std::vector<double>(n, INF));

    boost::johnson_all_pairs_shortest_paths(
        g, D,
        boost::weight_map(get(boost::edge_weight, g))
        .distance_inf(INF).distance_zero(0.0)
    );

    return D;
}

// ---------- bucket pairs by distance ----------
std::vector<std::vector<std::pair<int,int>>>
make_distance_buckets(const std::vector<std::vector<double>>& dist) {
    const int n = dist.size();
    double maxd = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (std::isfinite(dist[i][j]))
                maxd = std::max(maxd, dist[i][j]);

    int D = static_cast<int>(std::ceil(maxd));
    std::vector<std::vector<std::pair<int,int>>> buckets(D + 1);

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (std::isfinite(dist[i][j])) {
                int idx = static_cast<int>(std::round(dist[i][j]));
                buckets[std::min(idx, D)].emplace_back(i, j);
            }

    return buckets;
}

// ---------- sample N unique pairs ----------
std::vector<std::pair<int,int>> sample_unique_pairs_by_distance(
    const std::vector<std::vector<std::pair<int,int>>>& buckets,
    int N,
    std::mt19937& rng)
{
    std::vector<std::pair<int,int>> samples;
    samples.reserve(N);
    std::unordered_set<std::pair<int,int>, PairHash> seen;

    std::uniform_int_distribution<int> dsel(1, (int)buckets.size() - 1);

    while ((int)samples.size() < N) {
        int d = dsel(rng);
        if (buckets[d].empty()) continue;

        const auto& b = buckets[d];
        std::uniform_int_distribution<int> psel(0, (int)b.size() - 1);
        auto p = b[psel(rng)];

        if (seen.insert(p).second)
            samples.push_back(p);
    }

    return samples;
}



int main() {

    /*
    
    READING CSV FILE: EDGE[0], EDGE[1], WEIGHT

    */

    Graph read_g(0);

    std::ifstream file("/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/KSP/graphs_qspbwss/scale_free/BA_m2_n500.csv");
    std::string line;
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

    int K_max = 5000;

    auto dist = all_pairs_shortest_paths(read_g);
    auto buckets = make_distance_buckets(dist);

    std::mt19937 rng(std::random_device{}());
    int N = 300;
    auto samples = sample_unique_pairs_by_distance(buckets, N, rng);
    printf("Sampled %d unique pairs.\n", (int)samples.size());

    for (int i = 0; i < (int)samples.size(); ++i) {
        auto [s, t] = samples[i];

        std::vector<std::set<Vertex>> ini_p;
        std::vector<double> dist_ini;

        //double start = omp_get_wtime();

        std::pair<std::vector<Path>, std::vector<double>> pair_t = get_KSPs(read_g, s, t, ini_p, dist_ini, K_max);
        printf("Processed pair %d/%d: (%d, %d)\n", i, N, samples[i].first, samples[i].second);

        std::vector<Path> paths = std::get<0>(pair_t);
        std::vector<double> pathLens = std::get<1>(pair_t);


        //Writing file with paths and distance distribution
        std::ostringstream pathStream;
        // pathStream << "/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/KSP/graphs_qspbwss/osm_bcn_files/osm_bcn_distrs/path_distr/eix/paths_"
        //            << s << '_' << t << ".csv";
        pathStream << "/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/KSP/graphs_qspbwss/scale_free/paths_"
            << s << '_' << t << ".csv";
        std::ofstream fout(pathStream.str());

        std::ostringstream distStream;
        // distStream << "/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/KSP/graphs_qspbwss/osm_bcn_files/osm_bcn_distrs/path_distr/eix/dists_"
        //            << s << '_' << t << ".csv";
        distStream << "/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/KSP/graphs_qspbwss/scale_free/dists_"
            << s << '_' << t << ".csv";
        std::ofstream f2out(distStream.str());

        for (int j = 0; j < paths.size();++j) {
            const auto &path = paths[j];
            f2out << pathLens[j] << "\n";
            for (int i = 0; i < path.size(); ++i) {
                fout << path[i];
                if (i < path.size() - 1) {
                    fout << ",";
                }
            }
            fout << "\n";
        }
        fout.close();
        f2out.close();
    
    }
    return 0;

}
