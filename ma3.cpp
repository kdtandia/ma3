#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <limits>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>

using namespace std;

const int INF = numeric_limits<int>::max();

// Structure to represent an edge in the graph
struct Edge
{
    string destination;
    int weight;
};

// Graph class using adjacency list representation
class Graph
{
private:
    map<string, vector<Edge>> adjList;
    vector<string> cities;

public:
    // Add an edge to the graph
    void addEdge(const string &source, const string &dest, int weight)
    {
        adjList[source].push_back({dest, weight});

        // Keep track of all unique cities
        if (find(cities.begin(), cities.end(), source) == cities.end())
        {
            cities.push_back(source);
        }
        if (find(cities.begin(), cities.end(), dest) == cities.end())
        {
            cities.push_back(dest);
        }
    }

    // Get all cities in the graph
    vector<string> getCities() const
    {
        return cities;
    }

    // Get neighbors of a city
    const vector<Edge> &getNeighbors(const string &city) const
    {
        static const vector<Edge> empty;
        auto it = adjList.find(city);
        if (it != adjList.end())
        {
            return it->second;
        }
        return empty;
    }

    // Load graph from CSV file
    bool loadFromCSV(const string &filename)
    {
        ifstream file(filename);
        if (!file.is_open())
        {
            cerr << "Error: Could not open file " << filename << endl;
            return false;
        }

        string line;
        getline(file, line); // Skip header line

        while (getline(file, line))
        {
            stringstream ss(line);
            string source, dest, distStr;

            getline(ss, source, ',');
            getline(ss, dest, ',');
            getline(ss, distStr, ',');

            // Skip empty lines or lines with missing data
            if (source.empty() || dest.empty() || distStr.empty())
            {
                continue;
            }

            try
            {
                int distance = stoi(distStr);
                addEdge(source, dest, distance);
            }
            catch (...)
            {
                // Skip invalid lines
                continue;
            }
        }

        file.close();
        cout << "Loaded graph with " << cities.size() << " cities" << endl;
        return true;
    }

    // Brute Force Shortest Path using relaxation without priority queue
    // This repeatedly relaxes all edges until no more updates occur
    // Similar to Bellman-Ford but less efficient than Dijkstra
    map<string, int> bruteForceShortestPath(const string &start)
    {
        map<string, int> distances;

        // Initialize all distances to infinity
        for (const string &city : cities)
        {
            distances[city] = INF;
        }
        distances[start] = 0;

        // Relax edges repeatedly until no changes occur
        // This is inefficient because it doesn't use a priority queue
        bool changed = true;
        while (changed)
        {
            changed = false;

            // Try to relax all edges from all cities
            for (const string &city : cities)
            {
                if (distances[city] == INF)
                    continue;

                for (const Edge &edge : getNeighbors(city))
                {
                    int newDist = distances[city] + edge.weight;
                    if (newDist < distances[edge.destination])
                    {
                        distances[edge.destination] = newDist;
                        changed = true;
                    }
                }
            }
        }

        return distances;
    }

     map<string, int> dijkstraShortestPath(const string &start)
    {
        map<string, int> distances;

        // Initialize all distances to infinity
        for (const string &city : cities)
        {
            distances[city] = INF;
        }
        distances[start] = 0;

        
        priority_queue<pair<int, string>, vector<pair<int, string>>, greater<pair<int, string>>> pq;
        
        // Set to track visited nodes
        set<string> visited;
        
        // Start with the source city
        pq.push({0, start});
        
        while (!pq.empty())
        {
            // Get the city with the smallest distance
            string currentCity = pq.top().second;
            int currentDist = pq.top().first;
            pq.pop();
            
            // If we've already processed this city, skip it
            if (visited.find(currentCity) != visited.end())
            {
                continue;
            }
            
            // Mark as visited
            visited.insert(currentCity);
            
            // Check all neighbors of the current city
            for (const Edge &edge : getNeighbors(currentCity))
            {
                string neighbor = edge.destination;
                int weight = edge.weight;
                
                // If we haven't visited this neighbor and found a shorter path
                if (visited.find(neighbor) == visited.end())
                {
                    int newDist = currentDist + weight;
                    
                    // If we found a shorter path to the neighbor
                    if (newDist < distances[neighbor])
                    {
                        distances[neighbor] = newDist;
                        pq.push({newDist, neighbor});
                    }
                }
            }
        }

        return distances;
    }
};

// Print statistics (min, max, average, standard deviation)
void printStatistics(const string &algorithmName, const vector<double> &durations)
{
    double minTime = *min_element(durations.begin(), durations.end());
    double maxTime = *max_element(durations.begin(), durations.end());
    double avgTime = accumulate(durations.begin(), durations.end(), 0.0) / durations.size();

    double variance = 0.0;
    for (double t : durations)
    {
        variance += (t - avgTime) * (t - avgTime);
    }
    double stdDev = sqrt(variance / durations.size());

    cout << "\n=== " << algorithmName << " Statistics ===\n";
    cout << "Min Time: " << minTime << " ms\n";
    cout << "Max Time: " << maxTime << " ms\n";
    cout << "Average Time: " << avgTime << " ms\n";
    cout << "Standard Deviation: " << stdDev << " ms\n";
}

int main()
{
    // Load the graph from CSV
    Graph graph;
    if (!graph.loadFromCSV("state_capitals.csv"))
    {
        return 1;
    }

    vector<string> cities = graph.getCities();
    cout << "\nComputing shortest paths from each capital to all other capitals...\n";
    cout << "Total cities: " << cities.size() << endl;

    const int TRIALS = 10;

    // ========== BRUTE FORCE ALGORITHM BENCHMARKING ==========
    cout << "\n========== BRUTE FORCE ALGORITHM ==========\n";
    vector<double> bruteForceDurations;

    for (int trial = 0; trial < TRIALS; ++trial)
    {
        auto start = chrono::high_resolution_clock::now();

        // Compute shortest paths from each city to all other cities
        for (const string &city : cities)
        {
            map<string, int> distances = graph.bruteForceShortestPath(city);
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> elapsed = end - start;
        bruteForceDurations.push_back(elapsed.count());

        cout << "Trial " << trial + 1 << " duration: " << elapsed.count() << " ms\n";
    }

    printStatistics("Brute Force", bruteForceDurations);

    // ========== DIJKSTRA'S ALGORITHM BENCHMARKING ==========
    cout << "\n========== DIJKSTRA'S ALGORITHM ==========\n";
    vector<double> dijkstraDurations;

    for (int trial = 0; trial < TRIALS; ++trial)
    {
        auto start = chrono::high_resolution_clock::now();

        // Compute shortest paths from each city to all other cities
        for (const string &city : cities)
        {
            map<string, int> distances = graph.dijkstraShortestPath(city);
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> elapsed = end - start;
        dijkstraDurations.push_back(elapsed.count());

        cout << "Trial " << trial + 1 << " duration: " << elapsed.count() << " ms\n";
    }

    printStatistics("Dijkstra", dijkstraDurations);

    // Optional: Verify correctness by comparing a sample path
    cout << "\n========== SAMPLE SHORTEST PATH ==========\n";
    string sampleStart = "Olympia WA";
    cout << "Computing shortest paths from " << sampleStart << ":\n\n";

    map<string, int> bruteForceResult = graph.bruteForceShortestPath(sampleStart);

    // Display a few sample distances
    vector<string> sampleCities = {"Salem OR", "Boise ID", "Sacramento CA", "Phoenix AZ", "Austin TX"};
    cout << "Brute Force Results:\n";
    for (const string &city : sampleCities)
    {
        if (bruteForceResult.find(city) != bruteForceResult.end())
        {
            int dist = bruteForceResult[city];
            if (dist == INF)
            {
                cout << "  " << sampleStart << " -> " << city << ": No path\n";
            }
            else
            {
                cout << "  " << sampleStart << " -> " << city << ": " << dist << " miles\n";
            }
        }
    }

    // Compare with Dijkstra results
    map<string, int> dijkstraResult = graph.dijkstraShortestPath(sampleStart);
    cout << "\nDijkstra Results:\n";
    for (const string &city : sampleCities)
    {
        if (dijkstraResult.find(city) != dijkstraResult.end())
        {
            int dist = dijkstraResult[city];
            if (dist == INF)
            {
                cout << "  " << sampleStart << " -> " << city << ": No path\n";
            }
            else
            {
                cout << "  " << sampleStart << " -> " << city << ": " << dist << " miles\n";
            }
        }
    }

    // Verify that both algorithms produce the same results
    bool allMatch = true;
    for (const string &city : cities)
    {
        if (bruteForceResult[city] != dijkstraResult[city])
        {
            allMatch = false;
            cout << "\nMismatch found for " << city << "!\n";
            cout << "  Brute Force: " << bruteForceResult[city] << "\n";
            cout << "  Dijkstra: " << dijkstraResult[city] << "\n";
        }
    }

    if (allMatch)
    {
        cout << "\n✓ All shortest path distances match between both algorithms!\n";
    }

    return 0;
}