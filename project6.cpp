//----------------PROF-MATTA-CS340-SOPHIA-LICKO-PROJECT6------------------//
//------------------------------------------------------------------------//
//HEADERS//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>
#include <chrono>

//HEADERS//
using namespace std;
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//START OF CITY STRUCT//

struct City {
    double x, y;
};

//END OF CITY STRUCT//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//START OF STRUCT EDGE//

// Structure to represent an edge
struct Edge {
    int to;
    double weight;
};

//END OF STRUCT EDGE//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//START OF CLASS GRAPH//

class Graph {
public:
    map<int, vector<Edge>> adjList;

    void addNode(int vertex) {
        adjList[vertex];
    }

    void addEdge(int from, int to, double weight) {
        adjList[from].push_back({to, weight});
    }
};

//END OF CLASS GRAPH//
//------------------------------------------------------------------------//
//FUNC DEFINITIONS//

//function to calculate Euclidean distance between two cities
double calculateDistance(const City& a, const City& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

//function to calculate total distance of tour
double totalDistance(const vector<int>& tour, const vector<City>& cities) {
    double distance = 0;
    for (size_t i = 0; i < tour.size(); i++) {
        distance += calculateDistance(cities[tour[i]], cities[tour[(i + 1) % tour.size()]]);
    }
    return distance;
}

//function to read graph from file
vector<vector<Edge>> readGraphFromFile(const string& filename) {
  ifstream file(filename);
    string line;
    vector<vector<Edge>> graph;

    if (!file.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    while (getline(file, line)) {
        istringstream iss(line);
        int from, to;
        double weight;
        iss >> from; //read vertex num

        vector<Edge> edges;
        while (iss >> to >> weight) {
            edges.push_back({to, weight});
        }
        graph.push_back(edges);
    }

    file.close();
    return graph;
}

//function to read city coordinates from a csv file
vector<vector<City>> readCoordinatesForAllGraphs(const string& filename) {
    ifstream file(filename);
    string line;
    vector<vector<City>> allGraphsCoordinates(10); // 10 graphs

    if (!file.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    //skip the header line && xy line
    getline(file, line);
    getline(file, line);

    while (getline(file, line)) {
        istringstream iss(line);

        for (int i = 0; i < 10; ++i) { //read coordinates for each graph
            double x, y;
            iss >> x; //read x
            if (iss.peek() == ',') iss.ignore(); //skip comma
            iss >> y; //read y
            allGraphsCoordinates[i].push_back(City{x, y});
            if (i < 9) iss.ignore(numeric_limits<streamsize>::max(), ','); //skip to next pair
        }
    }

    file.close();
    return allGraphsCoordinates;
}

//function for brute forceeee
pair<vector<int>, double> bruteForceTSP(const vector<vector<Edge>>& graph, const vector<City>& cities) {

    //read check
       if (graph.size() != cities.size()) {
        cerr << "Error: Graph size does not match number of cities." << endl;
        exit(EXIT_FAILURE);
    }

    vector<int> bestTour, currentTour(cities.size());
    double minDistance = numeric_limits<double>::max();

    //populate currentTour with all cities
    iota(currentTour.begin(), currentTour.end(), 0);

    do {
        double currentDistance = 0;
        for (size_t i = 0; i < currentTour.size(); i++) {
            int from = currentTour[i], to = currentTour[(i + 1) % currentTour.size()];

            //calculate distance between cities
            currentDistance += calculateDistance(cities[from], cities[to]);
        }

        //update minDistance & bestTour if shorter tour found
        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            bestTour = currentTour;
        }
    } while (next_permutation(currentTour.begin(), currentTour.end()));

    return {bestTour, minDistance};
}

//Prim's alg for MST
vector<int> constructMST(const vector<vector<Edge>>& graph, const vector<City>& cities) {
    int n = cities.size();
    vector<bool> inMST(n, false);
    vector<double> minWeight(n, numeric_limits<double>::max());
    vector<int> parent(n, -1);
    set<pair<double, int>> minHeap;

    minWeight[0] = 0;
    minHeap.insert({0, 0});

    while (!minHeap.empty()) {
        int u = minHeap.begin()->second;
        minHeap.erase(minHeap.begin());
        inMST[u] = true;

        for (int v = 0; v < n; ++v) {
            double weight = calculateDistance(cities[u], cities[v]);
            if (!inMST[v] && weight < minWeight[v]) {
                minHeap.erase({minWeight[v], v});
                minWeight[v] = weight;
                minHeap.insert({minWeight[v], v});
                parent[v] = u;
            }
        }
    }

    return parent;
}

//convert parent array to adjacency list
vector<vector<int>> toAdjacencyList(const vector<int>& parent) {
    int n = parent.size();
    vector<vector<int>> adjList(n);
    for (int i = 1; i < n; ++i) {
        adjList[parent[i]].push_back(i);
        adjList[i].push_back(parent[i]);  //for undirected graph
    }
    return adjList;
}

//function for preorder traversal of MST
void preorderTraversal(int u, const vector<vector<int>>& adjList, vector<bool>& visited, vector<int>& tour) {
    visited[u] = true;
    tour.push_back(u);
    for (int v : adjList[u]) {
        if (!visited[v]) {
            preorderTraversal(v, adjList, visited, tour);
        }
    }
}

//function for preorder walk
vector<int> preorderWalk(const vector<vector<int>>& adjList) {
    vector<int> tour;
    vector<bool> visited(adjList.size(), false);
    preorderTraversal(0, adjList, visited, tour);
    return tour;
}

//MST-based TSP
pair<vector<int>, double> mstBasedTSP(const vector<vector<Edge>>& graph, const vector<City>& cities) {
    //update constructMST function to accept graph
    vector<int> mst = constructMST(graph, cities);
    vector<vector<int>> adjList = toAdjacencyList(mst);

    vector<int> tour = preorderWalk(adjList);
    double distance = 0;

    //calculate total distance using graph
    for (size_t i = 0; i < tour.size() - 1; i++) {
        int from = tour[i], to = tour[i + 1];
        auto it = find_if(graph[from].begin(), graph[from].end(), [to](const Edge& e) { return e.to == to; });
        if (it != graph[from].end()) {
            distance += it->weight;
        }
    }
    //add distance to return to start
    int last = tour.back(), start = tour.front();
    auto it = find_if(graph[last].begin(), graph[last].end(), [start](const Edge& e) { return e.to == start; });
    if (it != graph[last].end()) {
        distance += it->weight;
    }

    return {tour, distance};
}

//function for brute force alg to be timed
pair<vector<int>, double> timedBruteForceTSP(const vector<vector<Edge>>& graph, const vector<City>& cities) {
    auto startTime = chrono::high_resolution_clock::now();

    auto result = bruteForceTSP(graph, cities);

    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> execTime = endTime - startTime;
    cout << "Brute-Force TSP Execution Time: " << execTime.count() << " seconds." << endl;

    return result;
}

//function for MST alg to be timed
pair<vector<int>, double> timedMSTBasedTSP(const vector<vector<Edge>>& graph, const vector<City>& cities) {
    auto startTime = chrono::high_resolution_clock::now();

    auto result = mstBasedTSP(graph, cities);

    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> execTime = endTime - startTime;
    cout << "MST-Based TSP Execution Time: " << execTime.count() << " seconds." << endl;

    return result;
}


//FUNC DEFINITIONS//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//START OF MAIN//
int main() {
const string folderPath = "project6Graphs/";
    const string coordinatesFilename = folderPath + "coordinates.csv";

    //read coordinates for all graphs
    vector<vector<City>> allGraphsCoordinates = readCoordinatesForAllGraphs(coordinatesFilename);

    //iterate over all graphs
    for (int i = 0; i < allGraphsCoordinates.size(); ++i) {
        string graphFilename = folderPath + "graph" + to_string(i + 1) + ".txt";
        vector<vector<Edge>> graph = readGraphFromFile(graphFilename);

        //ensure graph & city coordinates sizes match
        if (graph.size() != allGraphsCoordinates[i].size()) {
            cerr << "Error: Graph size does not match number of cities for graph " << i + 1 << endl;
            continue; //skip to next graph
        }

        cout << "Processing Graph " << i + 1 << "..." << endl;

        //Brute Force TSP w/ timing
        auto bruteForceStart = chrono::high_resolution_clock::now();
        auto bruteForceResult = bruteForceTSP(graph, allGraphsCoordinates[i]);
        auto bruteForceEnd = chrono::high_resolution_clock::now();
        chrono::duration<double> bruteForceExecTime = bruteForceEnd - bruteForceStart;

        cout << "Graph " << i + 1 << " - Shortest tour by brute-force: ";
        for (int city : bruteForceResult.first) {
            cout << city << " ";
        }
        cout << "\nDistance (Brute Force): " << bruteForceResult.second << endl;
        cout << "Execution Time (Brute Force): " << bruteForceExecTime.count() << " seconds." << endl;

        //MST-based TSP w/ timing
        auto mstStart = chrono::high_resolution_clock::now();
        auto mstResult = mstBasedTSP(graph, allGraphsCoordinates[i]);
        auto mstEnd = chrono::high_resolution_clock::now();
        chrono::duration<double> mstExecTime = mstEnd - mstStart;

        cout << "Graph " << i + 1 << " - Approximate tour by MST: ";
        for (int city : mstResult.first) {
            cout << city << " ";
        }
        cout << "\nDistance (MST): " << mstResult.second << endl;
        cout << "Execution Time (MST): " << mstExecTime.count() << " seconds." << endl;
    }

    return 0;
}
//------------------------------------------------------------------------//
//----------------PROF-MATTA-CS340-SOPHIA-LICKO-PROJECT6------------------//