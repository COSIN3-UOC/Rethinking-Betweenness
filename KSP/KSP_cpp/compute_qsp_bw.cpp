#include "qsp_bw.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    int max_files = std::numeric_limits<int>::max();
    if (argc > 1) {
        try {
            max_files = std::stoi(argv[1]);
        } catch (const std::exception& ex) {
            std::cerr << "Invalid max-files argument: " << ex.what() << "\n";
            return 1;
        }
    }

    fs::path graph_folder = "../graphs_qspbwss/ring_roads/vx2/graph_variations";

    int file_count = 1;
    for (const auto& entry : fs::directory_iterator(graph_folder)) {
        std::cout << "Processing file: " << file_count << "/1000" << std::endl;
        if (entry.path().extension() != ".csv") {
            ++file_count;
            continue;
        }

        Graph read_g(0);
        std::ifstream file(entry.path());
        std::string line;

        if (!file.is_open()) {
            std::cerr << "Unable to open graph file: " << entry.path() << "\n";
            ++file_count;
            continue;
        }

        std::getline(file, line);  // skip header
        while (std::getline(file, line)) {
            boost::tokenizer<boost::escaped_list_separator<char>> tokens(line);
            auto tokenIterator = tokens.begin();
            Vertex vertex1 = std::stoi(*tokenIterator++);
            Vertex vertex2 = std::stoi(*tokenIterator++);
            double weight = std::stod(*tokenIterator);
            Edge e = add_edge(vertex1, vertex2, read_g).first;
            put(edge_weight, read_g, e, weight);
        }

        auto index = get(&VertexProperties::custom_index, read_g);
        graph_traits<Graph>::vertex_iterator vi, vend;
        graph_traits<Graph>::vertices_size_type cnt = 0;
        for (boost::tie(vi, vend) = vertices(read_g); vi != vend; ++vi) {
            put(index, *vi, cnt++);
        }

        std::vector<double> tol = {0.1, 0.2, 0.3};

        auto start = high_resolution_clock::now();
        auto betw = compute_tolSP_betw(read_g, tol);
        auto stop = high_resolution_clock::now();

        std::chrono::duration<float> duration = (stop - start);
        std::cout << "Time taken for this graph: "
                  << duration.count() << " seconds" << std::endl;

        std::string filename = entry.path().filename().string();
        int graph_number = extractFirstNumber(filename);
        std::string out_path = "../graphs_qspbwss/ring_roads/vx2_tst/bwss_" + std::to_string(graph_number);
        std::filesystem::create_directory(out_path);
        writeBTWtoCSV(betw, out_path, tol);

        ++file_count;
        if (file_count > max_files) {
            break;
        }
    }

    return 0;
}
