#pragma once
#include <vector>
#include <ostream>
#include <istream>
#include <random>

typedef std::vector<int> cluster_t;
typedef std::vector<cluster_t> clustering_t;
typedef std::vector<clustering_t> clusterings_t;

class Graph {
private:
	struct edge_struct {
		size_t node_start;
		size_t node_end;
		double cost;

		edge_struct(size_t node_start, size_t node_end, double cost) : node_start(node_start), node_end(node_end), cost(cost) {};
	};
	std::vector<std::string> V;
	std::vector<struct edge_struct> E;
	std::vector<std::vector<double>> E_matrix;
	size_t num_edges = 0;
	bool use_matrix = false;
public:
	// Calculate positive error for a given clustering
	// The positive error is the number of edges between nodes in different clusters
	double getPositiveErrorInClustering(const clustering_t& clusterings) const;

	// Calculate negative error for a given clustering
	// The negative error is the number of missing edges between nodes in the same cluster
	double getNegativeErrorInClustering(const clustering_t& clusterings) const;

	// Returns the sum of the positive and negative errors for the given clustering set
	double getErrorInClustering(const clustering_t& clusterings) const;

	// Streams a textual representation of the graph to the given output stream
	void streamGraph(std::ostream& stream) const;

	// Streams a textual representation of the graph to the given output stream, together with the clusterings
	void streamGraphClustering(std::ostream& stream, const clustering_t& clusterings) const;

	unsigned int getNumNodes() const;
	unsigned int getNumEdges() const;

	void addEdge(size_t index1, size_t index2);
	void addEdgeWithCost(size_t index1, size_t index2, double cost);
	void setNodeName(size_t index, std::string name);
	bool hasEdge(size_t index1, size_t index2) const;

	// Does not change existing edges to matrix, maybe should
	void useMatrixRepresentation(bool use);

	/**
		Creates a graph of the given size, with a "seed" (warning - not pseudorandom, but ordered)
		Warning: only works for graphs up to size approx. 11, afterwards will not generate all possible options
	**/
	Graph(int size, long long seed);

	/**
		Returns the number of possible different graphs of a certain size.
		Warning: only works for graphs up to size approx. 11 - otherwise overflows
	**/
	static long long numberOfPossibleGraphs(int size);

	static std::vector<Graph*> readG6Stream(std::istream& stream);
	static Graph* readJENAStream(std::istream& stream);

	// Gets the number of all possible edges in a graph of a particular size (no refexive edges)
	static int numberOfEdgesInAGraph(int size);

	/**
		Creates a graph of the given size, with a probability p of including each possible edge
	**/
	static Graph* createRandomGraph(int size, double p, std::mt19937_64& generator);

	/**
		Creates a graph of the given size, with a tendency to create clusters
			size - size of the graph to create
			cluster_size - size of the clusters
			cluster_size_stddev - standard deviation of the cluster_size variable (to generate clusters of different sizes)
			drop_rate - chance to drop edges in the generated graph
			add_rate - chance to add edges to the generated graph
	**/
	static Graph* createRandomClusterGraph(int size, double r, double r_stddev, double drop_rate, double add_rate, std::mt19937_64& generator);

	static Graph* createKPartiteCompleteGraph(std::vector<int> sizes);

	static clustering_t convertNodeListToClustering(int size, int* nodeToCluster);
};
