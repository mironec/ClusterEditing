#include "graph_utils.h"
#include <map>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>

double Graph::getPositiveErrorInClustering(const clustering_t& clusterings) const
{
	double positiveError = 0.0;

	std::map<int, int> nodeToCluster;

	int clusterNumber = 0;
	for (std::vector<int> cluster : clusterings) {
		for (int n : cluster) {
			nodeToCluster[n] = clusterNumber;
		}
		clusterNumber++;
	}

	if (use_matrix) {
		for (size_t i = 0; i < E_matrix.size(); i++) {
			for (size_t j = i+1; j < E_matrix[i].size(); j++) {
				assert(E_matrix.size() == E_matrix[i].size());
				if (nodeToCluster[i] != nodeToCluster[j] && E_matrix[i][j] > 0.0)					//Negative weights are assumed to be missing edges
					positiveError += E_matrix[i][j];
			}
		}
	}
	else {
		for (const struct edge_struct& edge : E) {
			if (nodeToCluster[edge.node_start] != nodeToCluster[edge.node_end] && edge.cost > 0.0)	//Negative weights are assumed to be missing edges
				positiveError += edge.cost;
		}
	}

	return positiveError;
}

double Graph::getNegativeErrorInClustering(const clustering_t& clusterings) const
{
	double negativeError = 0.0;

	for (const std::vector<int>& cluster : clusterings) {
		for (size_t index = 0; index < cluster.size(); index++) {
			for (size_t index2 = index + 1; index2 < cluster.size(); index2++) {		//Skip reflexive edges and do not double count
				int n1 = cluster[index];
				int n2 = cluster[index2];
				if (use_matrix) {
					if (E_matrix[n1][n2] < 0.0)
						negativeError -= E_matrix[n1][n2];
				}
				else {
					bool foundEdge = false;
					for (const struct edge_struct& edge : E) {
						if (((int)edge.node_start == n1 && (int)edge.node_end == n2) ||
							((int)edge.node_start == n2 && (int)edge.node_end == n1)) {			//Undirected graph assumed
							foundEdge = true;
							if (edge.cost < 0.0)										//Negative weights are assumed to be missing edges
								negativeError -= edge.cost;
							break;
						}
					}
					if (!foundEdge)
						negativeError++;
				}
			}
		}
	}

	return negativeError;
}

double Graph::getErrorInClustering(const clustering_t & clusterings) const
{
	return getPositiveErrorInClustering(clusterings) + getNegativeErrorInClustering(clusterings);
}

double Graph::getDeltaInClusteringMerge(const cluster_t & a, const cluster_t & b) const
{
	double deltaError = 0.0;

	if (use_matrix) {
		for (size_t n1 : a) {
			for (size_t n2 : b) {
				if (use_matrix) {
					// Works for all edges, negative, positive, 0-weight
					deltaError -= E_matrix[n1][n2];
				}
			}
		}
	}

	else { throw new std::exception("Not using matrix representation"); }

	return deltaError;
}

double Graph::getDeltaInClusteringSplit(const cluster_t & a, const cluster_t & b) const
{
	// Exactly inverse to the merge of the two clusters
	return -getDeltaInClusteringMerge(a, b);
}

void Graph::streamGraph(std::ostream & stream) const
{
	stream << "Graph" << std::endl;

	stream << "Nodes: ";
	for (std::string node_name : V) stream << node_name << ", ";
	stream << std::endl;

	stream << "Edges: ";
	if (use_matrix) {
		for (size_t i = 0; i < E_matrix.size(); i++) {
			for (size_t j = i+1; j < E_matrix[i].size(); j++) {
				if(E_matrix[i][j] != 0.0)
					stream << "(" << V[i] << ", " << V[j] << "): " << E_matrix[i][j] << ", ";
			}
		}
	}
	else {
		for (struct edge_struct e : E) stream << "(" << V[e.node_start] << ", " << V[e.node_end] << "): " << e.cost << ", ";
	}
	stream << std::endl;
}

void Graph::streamGraphClustering(std::ostream & stream, const clustering_t & clusterings) const
{
	stream << "Clusters { " << std::endl;
	int numCluster = 1;
	for (const std::vector<int>& cluster : clusterings) {
		stream << "  Cluster " << numCluster << " nodes: ";
		for (int n : cluster) stream << V[n] << ", ";
		stream << std::endl;
		numCluster++;
	}
	stream << "}" << std::endl;
}

unsigned int Graph::getNumNodes() const
{
	return V.size();
}

unsigned int Graph::getNumEdges() const
{
	if (use_matrix) return num_edges;
	return E.size();
}

void Graph::addEdge(size_t index1, size_t index2)
{
	addEdgeWithCost(index1, index2, 1.0);
}

void Graph::addEdgeWithCost(size_t index1, size_t index2, double cost)
{
	if (use_matrix) {
		if (E_matrix[index1][index2] != 0.0) num_edges--;
		E_matrix[index1][index2] = E_matrix[index2][index1] = cost;
		if(cost != 0.0)
			num_edges++;
	}
	else E.push_back(edge_struct(index1, index2, cost));
}

void Graph::setNodeName(size_t index, std::string name)
{
	V[index] = name;
}

bool Graph::hasEdge(size_t index1, size_t index2) const
{
	if (use_matrix) {
		return E_matrix[index1][index2] != 0.0;
	}
	else {
		for (const struct edge_struct& e : E) {
			if ((e.node_start == index1 && e.node_end == index2) || (e.node_start == index2 && e.node_end == index1)) return true;
		}
	}

	return false;
}

void Graph::useMatrixRepresentation(bool use)
{
	use_matrix = use;
	if (use) {
		assert(E_matrix.size() == 0);

		for (std::string v : V) {
			std::vector<double> vec;
			for (std::string v2 : V) {
				vec.push_back(0.0);
			}
			E_matrix.push_back(vec);
		}
	}
}

Graph::Graph(int size, long long seed)
{
	for (int i = 0; i < size; i++) {
		std::stringstream ss;
		ss << (i + 1);
		V.push_back(ss.str());
	}

	for (int n1 = 0; n1 < size; n1++) {
		for (int n2 = n1 + 1; n2 < size; n2++) {
			if (seed % 2 == 1) addEdge(n1, n2);
			seed = seed / 2;
		}
	}
}

long long Graph::numberOfPossibleGraphs(int size)
{
	return (long long)pow(2, size*(size - 1) / 2);
}

std::vector<Graph*> Graph::readG6Stream(std::istream& stream)
{
	std::vector<Graph*> graphs;

	int line = 0;
	int charnum = 0;
	int state = 0;
	long long current_graph_size = -1;
	int edge1 = 0; int edge2 = 1;
	bool ignore_rest_of_edges = false;
	Graph* current_graph = nullptr;

	while (stream.good()) {
		char c = stream.get();
		if (c == '>' && charnum == 0) {
			std::string header = ">>graph6<<";
			char readstring[10];
			stream.get(readstring, 9);
			if (header != readstring) {
				return graphs;
			}
			continue;
		}
		if (c == '\r') continue;
		if (c == '\n') {
			line++;
			state = 0;
			current_graph_size = -1;
			edge1 = 0;
			edge2 = 1;
			ignore_rest_of_edges = false;
			continue;
		}
		else if (state == 0) { // Read N(n) [one byte] - size of the graph
			if (c >= 63 && c <= 125) {
				current_graph_size = c - 63;
				state = 13;
			}
			else if (c == 126) {
				state = 1;
				current_graph_size = 0;
			}
			else {
				return graphs;
			}
		}
		else if (state >= 1 && state <= 4) { // Read N(n) [four bytes] - size of the graph
			if (c >= 63 && c <= 125) {
				current_graph_size += (c - 63) << ((4-state)*6);
				if (state < 4) state++;
				else state = 13;
			}
			else if (c == 126 && state == 1) {
				state = 5;
			}
			else {
				return graphs;
			}
		}
		else if (state >= 5 && state <= 12) { // Read N(n) [eight bytes] - size of the graph
			if (c >= 63 && c <= 125) {
				current_graph_size += (c - 63) << ((12-state)*6);
				if (state < 12) state++;
				else state = 13;
			}
			else {
				return graphs;
			}
		}

		else if (state == 13) {
			current_graph = new Graph((int)current_graph_size, 0);
			graphs.push_back(current_graph);

			state = 14;
		}

		if (state == 14) {
			if (c < 63) return graphs;
			if (ignore_rest_of_edges) continue;
			int edges = c - 63;

			for (int bitindex = 0; bitindex < 6; bitindex++) {
				int bit = (edges >> (5 - bitindex)) & 1;
				if (bit == 1) {
					current_graph->addEdge(edge1, edge2);
				}

				edge1++;
				if (edge1 == edge2) {
					edge2++;
					edge1 = 0;
					if (edge2 >= current_graph_size) {
						ignore_rest_of_edges = true;
						break;
					}
				}
			}
		}

		charnum++;
	}

	return graphs;
}

Graph * Graph::readJENAStream(std::istream & stream)
{
	std::string line;
	getline(stream, line);
	int graph_size = atoi(line.c_str());

	Graph* graph = new Graph(graph_size, 0);
	graph->useMatrixRepresentation(true);

	// Read node names
	for (int i = 0; i < graph_size; i++) {
		getline(stream, line);
		graph->setNodeName(i, line);
	}

	for (int i = 0; i < graph_size; i++) {
		getline(stream, line);
		size_t last_pos = -1;
		size_t pos = 0;
		int n1 = i;
		// Skip reflexive, and start from n1
		int n2 = i+1;
		while ((pos = line.find("\t", pos+1)) != std::string::npos) {
			std::string edge_desc = line.substr(last_pos+1, pos - last_pos - 1);
			double edge_cost = 0.0;

			// Bandaid fix for overflows
			if (edge_desc == "inf") edge_cost = 1e20;
			else if (edge_desc == "-inf") edge_cost = -1e20;
			else {
				edge_cost = atof(edge_desc.c_str());
				if (edge_cost < -1e250) {
					edge_cost = -1e20;
				}
				if (edge_cost > 1e250) {
					edge_cost = 1e20;
				}
			}

			graph->addEdgeWithCost(n1, n2, edge_cost);
			last_pos = pos;
			n2++;
		}
		std::string edge_desc = line.substr(last_pos + 1);
		if(edge_desc != "")
			graph->addEdgeWithCost(n1, n2, atof(edge_desc.c_str()));
	}

	return graph;
}

int Graph::numberOfEdgesInAGraph(int size)
{
	return size * (size - 1) / 2;
}

Graph * Graph::createRandomGraph(int size, double p, std::mt19937_64& generator)
{
	Graph* graph = new Graph(size, 0);

	int edge1 = 0;
	int edge2 = 1;
	int numEdges = numberOfEdgesInAGraph(size);
	for (long long i = 0; i < numEdges; i++) {
		if((double)generator() / generator.max() < p)
			graph->addEdge(edge1, edge2);
		edge2++;
		if (edge2 >= size) {
			edge1++;
			edge2 = edge1 + 1;
		}
	}

	return graph;
}

Graph * Graph::createRandomClusterGraph(int size, double cluster_size, double cluster_size_stddev, double drop_rate, double add_rate, std::mt19937_64 & generator)
{
	Graph* graph = new Graph(size, 0);

	std::vector<int> permutation;
	permutation.reserve(size);
	for (int i = 0; i < size; i++) permutation.insert(permutation.begin() + (permutation.size() > 0 ? generator() % permutation.size() : 0), i);

	std::normal_distribution<double> cluster_size_distribution(cluster_size, cluster_size_stddev);
	std::vector<int> cluster_sizes;
	for (int used = 0; used < size; ) {
		int current_size = (int)cluster_size_distribution(generator);
		if (current_size <= 0) continue;
		if (current_size > size - used) current_size = size - used;

		cluster_sizes.push_back(current_size);
		used += current_size;
	}

	for (int n1 = 0; n1 < size; n1++) {
		for (int n2 = n1 + 1; n2 < size; n2++) {
			if (generator() / (double)generator.max() < drop_rate) continue;

			int n1_real = permutation[n1];
			int n2_real = permutation[n2]; if (n2_real < n1_real) std::swap(n1_real, n2_real);

			if (generator() / (double)generator.max() < add_rate) {
				graph->addEdge(n1_real, n2_real);
				continue;
			}

			int used = 0;
			for (int csize : cluster_sizes) {
				if (n1 < used + csize) {
					if (n2 >= used && n2 < used + csize) {
						graph->addEdge(n1_real, n2_real);
					}
					break;
				}
				used += csize;
			}
		}
	}

	return graph;
}

Graph * Graph::createKPartiteCompleteGraph(std::vector<int> sizes)
{
	int size_total = 0;
	std::vector<int> cumulative_sizes;
	for (int s : sizes) {
		cumulative_sizes.push_back(size_total);
		size_total += s;
	}
	Graph* graph = new Graph(size_total, 0);

	for (size_t i = 0; i < sizes.size(); i++) {
		for (size_t j = i+1; j < sizes.size(); j++) {
			for (int n1 = cumulative_sizes[i]; n1 < sizes[i] + cumulative_sizes[i]; n1++) {
				for (int n2 = cumulative_sizes[j]; n2 < sizes[j] + cumulative_sizes[j]; n2++) {
					graph->addEdge(n1, n2);
				}
			}
		}
	}

	return graph;
}

clustering_t Graph::convertNodeListToClustering(int size, int* nodeToCluster)
{
	clustering_t clustering;

	std::map<int, int> clusterNumberToIndex;
	for (int i = 0; i < size; i++) {
		int clusterNumber = nodeToCluster[i];
		if (clusterNumberToIndex.count(clusterNumber) > 0) {
			int index = clusterNumberToIndex[clusterNumber];
			clustering[index].push_back(i);
		}
		else {
			clusterNumberToIndex[clusterNumber] = clustering.size();
			std::vector<int> newvec;
			newvec.push_back(i);
			clustering.push_back(newvec);
		}
	}

	return clustering;
}
