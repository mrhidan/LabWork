#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string> 
#include <fstream>

class Graph {
public:
    using Vertex = size_t;
    using VertexSet = std::unordered_set<Vertex>;
    using AdjencyList = std::unordered_map<Vertex, VertexSet>;

    void AddVertex(Vertex v) {
        adjency_list_[v];
    }

    void AddEdge(Vertex u, Vertex v) {
        adjency_list_[u].insert(v);
        adjency_list_[v].insert(u);
    }

    void RemoveVertex(Vertex v) {
        adjency_list_.erase(v);
    }

    void RemoveEdge(Vertex u, Vertex v) {
        adjency_list_[u].erase(v);
        adjency_list_[v].erase(u);
    }

    const VertexSet& AdjecentVertices(Vertex v) const {
        const auto it = adjency_list_.find(v);
        if (it != adjency_list_.end()) {
            return it->second;
        }
        else {
            return empty_set_;
        }
    }

    VertexSet AllVertices() const {
        VertexSet vs;
        vs.reserve(adjency_list_.size());
        for (const auto& pair : adjency_list_) {
            const auto& vertex = pair.first;
            vs.insert(vertex);
        }
        return vs;
    }

    const AdjencyList& AsAdjencyList() const {
        return adjency_list_;
    }

private:
    AdjencyList adjency_list_;
    static const VertexSet empty_set_;
};

const Graph::VertexSet Graph::empty_set_;

class MaxLengthPath {
public:
    explicit  MaxLengthPath(const Graph& graph, Graph::Vertex startingpoint)
        : graph_(graph), set_({ startingpoint }), starting_point(startingpoint),
        FirstFinishNodeEncountered(graph.AdjecentVertices(startingpoint).size() > 2),
        FirstStartNodeEncountered(graph.AdjecentVertices(startingpoint).size() > 2) {}

    Graph::VertexSet CandidatesToAdd() const {
        Graph::VertexSet candidates;
        for (const auto v : graph_.AdjecentVertices(*path_start)) {
            if (IsAddable(v)) {
                candidates.insert(v);
            }
        }
        for (const auto v : graph_.AdjecentVertices(*path_finish)) {
            if (IsAddable(v)) {
                candidates.insert(v);
            }
        }
        return candidates;
    }

    void Add(Graph::Vertex v) {
        waygraph_.AddVertex(v);
        set_.insert(v);
        for (const auto u : graph_.AdjecentVertices(v)) {
            if (find(u) == finish()) {
                if (graph_.AdjecentVertices(v).size() > 2) {
                    if (!FirstFinishNodeEncountered) {
                        FirstFinishNodeDistance = 0;
                        FirstFinishNodeEncountered = true;
                    }
                    LastFinishNodeDistance = 0;
                }
                else {
                    ++LastFinishNodeDistance;
                }
                ++FirstFinishNodeDistance;
                waygraph_.AddVertex(v);
                waygraph_.AddEdge(v, *path_finish);
                path_finish = find(v);
                break;
            }
            if (find(u) == start()) {
                if (graph_.AdjecentVertices(v).size() > 2) {
                    if (!FirstStartNodeEncountered) {
                        FirstStartNodeDistance = 0;
                        FirstStartNodeEncountered = true;
                    }
                    LastStartNodeDistance = 0;
                }
                else {
                    ++LastStartNodeDistance;
                }
                ++FirstStartNodeDistance;
                waygraph_.AddVertex(v);
                waygraph_.AddEdge(v, *path_start);
                path_start = find(v);
                break;
            }
        }
    }

    void ReverseStart() {
        if (!waygraph_.AdjecentVertices(*path_start).empty()) {
            Graph::Vertex parent = *waygraph_.AdjecentVertices(*path_start).begin();
            if (*start() != starting_point) {
                size_t tmp = *path_start;
                set_.erase(tmp);
                waygraph_.RemoveEdge(tmp, *waygraph_.AdjecentVertices(tmp).begin());
                waygraph_.RemoveVertex(tmp);
                path_start = find(parent);
                LastStartNodeDistance -= LastStartNodeDistance > 0;
                --FirstStartNodeDistance;
            }
        }
    }

    void ReverseFinish() {
        if (!waygraph_.AdjecentVertices(*path_finish).empty()) {
            Graph::Vertex parent = *waygraph_.AdjecentVertices(*path_finish).begin();
            if (*finish() != starting_point) {
                size_t tmp = *path_finish;
                set_.erase(tmp);
                waygraph_.RemoveEdge(tmp, *waygraph_.AdjecentVertices(tmp).begin());
                waygraph_.RemoveVertex(tmp);
                path_finish = find(parent);
                LastFinishNodeDistance -= LastFinishNodeDistance > 0;
                --FirstFinishNodeDistance;
            }
        }
    }

    size_t length() const {
        return set_.size();
    }

    void operator= (MaxLengthPath& other) {
        waygraph_ = other.waygraph_;
        set_ = other.set_;
        path_start = other.path_start;
        path_finish = other.path_finish;
        FirstStartNodeDistance = other.FirstStartNodeDistance;
        FirstFinishNodeDistance = other.FirstFinishNodeDistance;
        LastFinishNodeDistance = other.LastFinishNodeDistance;
        LastStartNodeDistance = other.LastStartNodeDistance;
        FirstFinishNodeEncountered = other.FirstFinishNodeEncountered;
        FirstStartNodeEncountered = other.FirstStartNodeEncountered;
    }

    Graph::VertexSet::const_iterator find(Graph::Vertex v) const {
        return set_.find(v);
    }

    Graph::VertexSet::const_iterator begin() const {
        return set_.begin();
    }

    Graph::VertexSet::const_iterator end() const {
        return set_.end();
    }

    Graph::VertexSet::const_iterator start() const {
        return path_start;
    }

    Graph::VertexSet::const_iterator finish() const {
        return path_finish;
    }

    size_t GetLastStartNodeDistance() const {
        return LastStartNodeDistance;
    }

    size_t GetLastFinishNodeDistance() const {
        return LastFinishNodeDistance;
    }

    size_t GetFirstStartNodeDistance() const {
        return FirstStartNodeDistance;
    }

    size_t GetFirstFinishNodeDistance() const {
        return FirstFinishNodeDistance;
    }

    const Graph& GetGraph() const {
        return waygraph_;
    }

    bool IsFirstStartNodeEncountered() const {
        return FirstStartNodeEncountered;
    }

    bool IsFirstFinishNodeEncountered() const {
        return FirstFinishNodeEncountered;
    }

private:
    bool IsAddable(Graph::Vertex v) const {
        if (find(v) != end()) {
            return false;
        }
        return true;
    }
    Graph waygraph_;
    const Graph& graph_;
    Graph::VertexSet set_;
    Graph::Vertex starting_point;
    Graph::VertexSet::iterator path_start = set_.begin();
    Graph::VertexSet::iterator path_finish = set_.begin();
    size_t FirstStartNodeDistance = 0;
    size_t FirstFinishNodeDistance = 0;
    size_t LastStartNodeDistance = 0;
    size_t LastFinishNodeDistance = 0;
    bool FirstStartNodeEncountered;
    bool FirstFinishNodeEncountered;
};

void GraphEdges(std::ostream& out, const Graph::AdjencyList& adjency_list) {
    for (const auto& pair : adjency_list) {
        const auto& vertex = pair.first;
        const auto& neighbours = pair.second;
        for (const auto adj_vertex : neighbours) {
            out << "\t" << vertex << " -- " << adj_vertex << "\n";
        }
    }
}

// Use http://www.webgraphviz.com to take a look at the graph
void GraphViz(std::ostream& out, const Graph& graph) {
    out << "strict graph {\n";
    for (const auto& pair : graph.AsAdjencyList()) {
        const auto& vertex = pair.first;
        out << "\t" << vertex << "\n";
    }
    GraphEdges(out, graph.AsAdjencyList());
    out << "}\n";
}

void GraphViz(std::ostream& out, const MaxLengthPath& Max_Length_Path) {
    out << "strict graph {\n";
    for (const auto& pair : Max_Length_Path.GetGraph().AsAdjencyList()) {
        const auto& vertex = pair.first;
        if (Max_Length_Path.find(vertex) != Max_Length_Path.end()) {
            out << "\t" << vertex << " [shape=doublecircle]\n";
        }
        else {
            out << "\t" << vertex << "\n";
        }
    }
    GraphEdges(out, Max_Length_Path.GetGraph().AsAdjencyList());
    out << "}\n";
}

struct DebugInfo {
    std::vector<size_t> costs;
    std::size_t offset = 0;
};

// Use http://gnuplot.respawned.com/ to plot costs
std::ostream& operator<<(std::ostream& out, const DebugInfo& debug_info) {
    for (size_t i = 0; i < debug_info.costs.size(); ++i) {
        out << "\n" << i + debug_info.offset << " " << debug_info.costs[i];
    }
    return out;
}

class MaxLengthPathSolver {
public:
    virtual MaxLengthPath Solve(const Graph& graph,
        DebugInfo& debug_info) const = 0;
    virtual ~MaxLengthPathSolver() = default;
};

class GradientDescent final : public MaxLengthPathSolver { // Standard implementation
    MaxLengthPath Solve(const Graph& graph, DebugInfo& debug_info) const {
        size_t startingpoint = rand() % graph.AsAdjencyList().size();
        MaxLengthPath mlp = MaxLengthPath(graph, startingpoint);
        while (true) {
            debug_info.costs.push_back(mlp.length());
            auto cands = mlp.CandidatesToAdd();
            auto cand_it = cands.begin();
            if (cand_it == cands.end()) {
                break;
            }
            if (cands.size() == 1) {
                mlp.Add(*cand_it);
            }
            else {
                for (size_t i = 0; i < rand() % cands.size(); ++i) {
                    ++cand_it;
                }
                mlp.Add(*cand_it);
            }
        }
        return mlp;
    }
};

class PseudoMetropolis final : public MaxLengthPathSolver { // Implementation based on idea of optimizing backtrack vertex number
private:
    double k, T, cooling_rate;
    bool annealing;
public:
    PseudoMetropolis(double k_constant, double Temperature, bool IsAnnealing, double CoolingRate) {
        k = k_constant, T = Temperature, annealing = IsAnnealing,
            cooling_rate = 1 - CoolingRate; // In percentage
    }
    MaxLengthPath Solve(const Graph& graph, DebugInfo& debug_info) const {
        size_t startingpoint = rand() % graph.AsAdjencyList().size();
        auto operating_T = T;
        MaxLengthPath mlp = MaxLengthPath(graph, startingpoint);
        MaxLengthPath best_path = mlp;
        bool comparing = false;
        while (T != 0) {
            debug_info.costs.push_back(mlp.length());
            auto cands = mlp.CandidatesToAdd();
            auto cand_it = cands.begin();
            if (cand_it == cands.end()) {
                best_path = mlp;
                if (mlp.IsFirstFinishNodeEncountered() || mlp.IsFirstStartNodeEncountered()) {
                    bool reversing_start = rand() % static_cast <size_t> (2) == 1;
                    size_t r; // The number of vertex to give up
                    if ((mlp.IsFirstStartNodeEncountered()) && (reversing_start || (!mlp.IsFirstFinishNodeEncountered()))) {
                        r = mlp.GetFirstStartNodeDistance() != mlp.GetLastStartNodeDistance() ?
                            ((rand() + rand()) % (mlp.GetFirstStartNodeDistance() - mlp.GetLastStartNodeDistance()) + mlp.GetLastStartNodeDistance() + 1) : mlp.GetLastStartNodeDistance();
                        for (int i = 0; i < r; ++i) {
                            mlp.ReverseStart();
                        }
                    }
                    else if (mlp.IsFirstFinishNodeEncountered()) {
                        r = mlp.GetFirstFinishNodeDistance() != mlp.GetLastFinishNodeDistance() ?
                            ((rand() + rand()) % (mlp.GetFirstFinishNodeDistance() - mlp.GetLastFinishNodeDistance()) + mlp.GetLastFinishNodeDistance() + 1) : mlp.GetLastFinishNodeDistance();
                        for (int i = 0; i < r; ++i) {
                            mlp.ReverseFinish();
                        }
                    }
                }
                else {
                    return best_path;
                }
                double raw_1 = double(rand()) / (RAND_MAX); //Debug info
                double tmp = (static_cast <int> (mlp.length()) - static_cast <int> (best_path.length())) / static_cast <double> (k*operating_T); //Debug info
                double raw_2 = exp(tmp); //Debug info
                if (raw_1 > raw_2) {
                    return best_path;
                }
                else if (annealing) {
                    operating_T *= cooling_rate;
                }
            }
            else if (cands.size() == 1) {
                mlp.Add(*cand_it);
            }
            else {
                for (size_t i = 0; i < rand() % cands.size(); ++i) {
                    ++cand_it;
                }
                mlp.Add(*cand_it);
            }
            if (mlp.length() > best_path.length()) {
                best_path = mlp;
            }
        }
        return mlp;
    }
};

class Metropolis final : public MaxLengthPathSolver { // Implementation based on pure lecture material
private:
    double k, T, cooling_rate;
    bool annealing;
public:
    Metropolis(double k_constant, double Temperature, bool IsAnnealing, double CoolingRate) {
        k = k_constant, T = Temperature, annealing = IsAnnealing,
            cooling_rate = 1 - CoolingRate; // In percentage
    }
    MaxLengthPath Solve(const Graph& graph, DebugInfo& debug_info) const {
        size_t startingpoint = rand() % graph.AsAdjencyList().size();
        auto operating_T = T;
        MaxLengthPath mlp = MaxLengthPath(graph, startingpoint);
        bool comparing = false;
        while (T != 0) {
            debug_info.costs.push_back(mlp.length());
            auto cands = mlp.CandidatesToAdd();
            auto cand_it = cands.begin();
            if (cand_it == cands.end()) {
                double raw_1 = double(rand()) / (RAND_MAX); //Debug info
                double tmp = (-1) / static_cast <double> (k*operating_T); //Debug info
                double raw_2 = exp(tmp); //Debug info
                if (raw_1 < raw_2) {
                    if (rand() % static_cast <size_t> (2) == 1 && (mlp.IsFirstStartNodeEncountered())) {
                        mlp.ReverseStart();
                    }
                    else if (mlp.IsFirstFinishNodeEncountered()) {
                        mlp.ReverseFinish();
                    }
                }
                else {
                    break;
                }
            }
            else {
                size_t tmp = rand() % (cands.size() + 2);
                if (tmp < cands.size()) { // We want to add new vertex
                    for (size_t i = 0; i < tmp; ++i) {
                        ++cand_it;
                    }
                    mlp.Add(*cand_it);
                }
                else { // We want to backtrack
                    double raw_1 = double(rand()) / (RAND_MAX); //Debug info
                    double tmp = (-1) / static_cast <double> (k*operating_T); //Debug info, the difference between values is constant
                    double raw_2 = exp(tmp); //Debug info
                    if (raw_1 < raw_2) {
                        if (rand() % static_cast <size_t> (2) == 1 && (mlp.IsFirstStartNodeEncountered())) {
                            mlp.ReverseStart();
                        }
                        else if (mlp.IsFirstFinishNodeEncountered()) {
                            mlp.ReverseFinish();
                        }
                    }
                }
            }
            if (annealing) {
                operating_T *= cooling_rate;
            }
        }
        return mlp;
    }
};

Graph RandomGraph(size_t size, double edge_probability) {
    Graph graph;
    for (Graph::Vertex v = 1; v <= size; ++v) {
        graph.AddVertex(v);
    }
    for (Graph::Vertex v = 1; v <= size; ++v) {
        for (Graph::Vertex u = v + 1; u <= size; ++u) {
            if (double(rand() + rand()) / RAND_MAX <= edge_probability) {
                graph.AddEdge(v, u);
            }
        }
    }
    return graph;
}

Graph StarGraph(size_t size) {
    Graph graph;
    for (Graph::Vertex v = 2; v <= size; ++v) {
        graph.AddEdge(1, v);
    }
    return graph;
}

int InitRandSeed(int argc, const char* argv[]) {
    int rand_seed;
    if (argc >= 2) {
        rand_seed = atoi(argv[1]);
    }
    else {
        rand_seed = time(nullptr);
    }
    srand(rand_seed);
    return rand_seed;
}

void TrySolver(const MaxLengthPathSolver& solver, const Graph& graph, std::string methodname) {
    //GraphViz(std::cout, graph);
    std::ofstream datafile("data" + methodname + ".txt");
    datafile << "# Iterations Length Results";
    auto best_length = 0;
    size_t offset = 0;
    Graph best_length_graph;
    for (int attempt = 1; attempt < 100; ++attempt) {
        DebugInfo debug_info;
        const auto Max_Length_Path = solver.Solve(graph, debug_info);
        auto length = Max_Length_Path.length();
        debug_info.offset = offset;
        datafile << debug_info;
        offset = debug_info.offset + debug_info.costs.size();
        if (length > best_length) {
            best_length_graph = Max_Length_Path.GetGraph();
            best_length = length;
            //GraphViz(std::cout, Max_Length_Path);
            datafile << " " << best_length;
        }
    }
    datafile.close();
    //GraphViz(std::cout, best_length_graph);
    std::cout << "Result: " << best_length << std::endl;
}

int main(int argc, const char* argv[]) {
    std::cout << "Using rand seed: " << InitRandSeed(argc, argv) << "\n";
    const auto graph = RandomGraph(100, 0.21);
    std::ofstream graphdata("graph.txt");
    //const auto graph = StarGraph(100);
    GraphViz(graphdata, graph);
    graphdata.close();
    GradientDescent gradient_descent;
    PseudoMetropolis metropolis_1(1, 10, false, 0.05);
    PseudoMetropolis metropolis_2(1, 1000, true, 0.01);
    Metropolis metropolis_3(1, 10, false, 0.01);
    Metropolis metropolis_4(1, 1000, true, 0.01);
    TrySolver(metropolis_1, graph, "1.1");
    TrySolver(metropolis_2, graph, "1.2");
    TrySolver(gradient_descent, graph, "2");
    TrySolver(metropolis_3, graph, "3.1");
    TrySolver(metropolis_4, graph, "3.2");
    return 0;
}