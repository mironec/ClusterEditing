#include <iostream>
#include <fstream>
#include <sstream>
#include <experimental/filesystem>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <random>
#include <chrono>
#include <cassert>
#include <climits>
#include "graph_utils.h"
#include "Evaluator.h"
#include "GeneticAlgo.h"

using namespace std;

vector<vector<vector<int>>> bruteForce_bestClusterings;
vector<vector<int>> bruteForce_internal(const Graph& graph, int k, bool keepClusterings = false) {
	assert(k > 0);
	if (keepClusterings) bruteForce_bestClusterings.clear();

	map<int, int> nodeIndexToCluster;
	for (size_t i = 0; i < graph.getNumNodes(); i++) nodeIndexToCluster[i] = 0;

	bool moreClusters = true;
	double minError = INFINITY;
	vector<vector<int>> bestClusterings;
	while (moreClusters) {
		//Disqualify empty clusterings
		bool emptyClustering = false;
		int* numNodesInCluster = new int[k];
		for (int i = 0; i < k; i++) numNodesInCluster[i] = 0;
		for (pair<int, int> p : nodeIndexToCluster) {
			numNodesInCluster[p.second]++;
		}
		for (int i = 0; i < k; i++) if (numNodesInCluster[i] == 0) emptyClustering = true;

		if (!emptyClustering) {
			vector<vector<int>> clusterings;
			for (int i = 0; i < k; i++) clusterings.push_back(vector<int>());
			for (pair<int, int> p : nodeIndexToCluster) {
				clusterings[p.second].push_back(p.first);
			}

			double error = graph.getErrorInClustering(clusterings);
			if (error < minError) {
				bestClusterings = clusterings;
				minError = error;
				if (keepClusterings) bruteForce_bestClusterings.clear();
				if (keepClusterings) bruteForce_bestClusterings.push_back(clusterings);
			}
			else if (error == minError && keepClusterings) {
				bruteForce_bestClusterings.push_back(clusterings);
			}
		}

		//Make new configuration of clusterings
		moreClusters = false;
		for (size_t i = 0; i < graph.getNumNodes(); i++) {
			nodeIndexToCluster[i]++;
			if (nodeIndexToCluster[i] >= k) nodeIndexToCluster[i] = 0;
			else {
				moreClusters = true;
				break;
			}
		}

		delete[] numNodesInCluster;
	}

	return bestClusterings;
}

// k - number of clusterings to produce
vector<vector<int>> bruteForce(const Graph& graph, int k) {
	return bruteForce_internal(graph, k, false);
}

// k - number of clusterings to produce
vector<vector<int>> bruteForceKeepClusters(const Graph& graph, int k) {
	return bruteForce_internal(graph, k, true);
}

vector<vector<int>> bruteForce_internal_all(const Graph& graph) {
	unsigned int graph_size = graph.getNumNodes();

	map<unsigned int, unsigned int> nodeIndexToCluster;
	for (size_t i = 0; i < graph_size; i++) nodeIndexToCluster[i] = 0;

	bool moreClusters = true;
	double minError = INFINITY;
	vector<vector<int>> bestClusterings;
	while (moreClusters) {
		vector<vector<int>> clusterings;
		for (unsigned int i = 0; i < graph_size; i++) clusterings.push_back(vector<int>());
		for (pair<int, int> p : nodeIndexToCluster) {
			clusterings[p.second].push_back(p.first);
		}

		double error = graph.getErrorInClustering(clusterings);
		if (error < minError) {
			bestClusterings = clusterings;
			minError = error;
		}

		//Make new configuration of clusterings
		moreClusters = false;
		for (size_t i = 0; i < graph.getNumNodes(); i++) {
			nodeIndexToCluster[i]++;
			if (nodeIndexToCluster[i] >= graph_size) nodeIndexToCluster[i] = 0;
			else {
				moreClusters = true;
				break;
			}
		}
	}

	return bestClusterings;
}

enum EDGE_STATUS : char
{
	INITIALLY_PRESENT = 1,
	INITIALLY_MISSING = 2,
	PRESENT = 4,
	MISSING = 8,
	PERMANENT = 16,
	FORBIDDEN = 32
};

inline EDGE_STATUS operator|(EDGE_STATUS a, EDGE_STATUS b) { return static_cast<EDGE_STATUS>(static_cast<int>(a) | static_cast<int>(b)); }
inline EDGE_STATUS operator&(EDGE_STATUS a, EDGE_STATUS b) { return static_cast<EDGE_STATUS>(static_cast<int>(a) & static_cast<int>(b)); }

clustering_t findOptimalClusteringConflictTriples(const Graph& graph, int k) {
	assert(k > 0);

	int graph_size = graph.getNumNodes();
	int edge_size = Graph::numberOfEdgesInAGraph(graph_size);

	struct modification_node {
		modification_node(EDGE_STATUS* edge_status, int cost) : edge_status(edge_status), cost(cost) {};
		EDGE_STATUS *edge_status;
		int cost;
		struct modification_node* left = nullptr;				//complete the triangle
		struct modification_node* middle = nullptr;				//delete one branch, set the other to permanent, mark rest of triangle forbidden
		struct modification_node* right = nullptr;				//delete the other branch, set just that edge to forbidden
	};

	EDGE_STATUS *initial_edge_status = new EDGE_STATUS[edge_size];
	int* n1s = new int[edge_size];
	int* n2s = new int[edge_size];
	int** n1_n2_index = new int*[graph_size];
	for (int i = 0; i < graph_size; i++) n1_n2_index[i] = new int[graph_size];

	size_t index = 0;
	for (int n1 = 0; n1 < graph_size; n1++) {
		for (int n2 = n1 + 1; n2 < graph_size; n2++) {
			if (graph.hasEdge(n1, n2)) initial_edge_status[index] = EDGE_STATUS::INITIALLY_PRESENT | EDGE_STATUS::PRESENT;
			else initial_edge_status[index] = EDGE_STATUS::INITIALLY_MISSING | EDGE_STATUS::MISSING;
			n1s[index] = n1;
			n2s[index] = n2;
			n1_n2_index[n1][n2] = index;
			index++;
		}
	}
	modification_node root(initial_edge_status, 0);

	EDGE_STATUS *current_edge_status = initial_edge_status;

	for (int n1 = 0; n1 < graph_size; n1++) {
		for (int n2 = n1 + 1; n2 < graph_size; n2++) {
			for (int n3 = n2 + 1; n3 < graph_size; n3++) {
				// Bad triplet
				if ((current_edge_status[n1_n2_index[n1][n2]] & EDGE_STATUS::PRESENT) &&
					(current_edge_status[n1_n2_index[n1][n3]] & EDGE_STATUS::PRESENT) &&
					(current_edge_status[n1_n2_index[n2][n3]] & EDGE_STATUS::MISSING) &&
					(!(current_edge_status[n1_n2_index[n1][n2]] & EDGE_STATUS::PERMANENT) ||
					!(current_edge_status[n1_n2_index[n1][n3]] & EDGE_STATUS::PERMANENT) ||
					!(current_edge_status[n1_n2_index[n2][n3]] & EDGE_STATUS::FORBIDDEN)) ) {
					//pass
				}
			}
		}
	}

	delete[] initial_edge_status;
	delete[] n1s;
	delete[] n2s;

	return {};
}

bool isNegativeBond(const Graph& graph, const cluster_t& cluster1, const cluster_t& cluster2) {
	size_t c1_size = cluster1.size();
	size_t c2_size = cluster2.size();
	size_t comb_size = c1_size + c2_size;
	int bond_strength = Graph::numberOfEdgesInAGraph(c1_size) + Graph::numberOfEdgesInAGraph(c2_size) - Graph::numberOfEdgesInAGraph(comb_size);

	for (int n1 : cluster1) {
		for (int n2 : cluster2) {
			if (graph.hasEdge(n1, n2)) {
				bond_strength += 2; // The above bond strength should be divided by 2, but I don't want fractions, so instead we multiply by 2 here
			}
		}
	}

	return bond_strength < 0;
}

bool hasAllNegativeBonds(const Graph& graph, const clustering_t& clustering) {
	size_t num_clusters = clustering.size();

	for (size_t c1 = 0; c1 < num_clusters; c1++) {
		for (size_t c2 = c1 + 1; c2 < num_clusters; c2++) {
			if (!isNegativeBond(graph, clustering[c1], clustering[c2])) return false;
		}
	}

	return true;
}

void findClustering(const Graph& graph, bool print_mode = false) {
	int size = graph.getNumNodes();
	double lastError = -1;
	//double lastPositiveError = -1;
	//double lastNegativeError = -1;
	bool achievedMinimum = false;

	for (int i = 0; i < size; i++) {
		vector<vector<int>> bestClustering = bruteForce(graph, i + 1);

		double pError = graph.getPositiveErrorInClustering(bestClustering);
		double nError = graph.getNegativeErrorInClustering(bestClustering);
		double error = pError + nError;
		/*for (const clustering_t& clustering : bruteForce_bestClusterings) {
			double p = graph.getPositiveErrorInClustering(clustering);
			double n = graph.getNegativeErrorInClustering(clustering);

			if (n < nError) { pError = p; nError = n; }
		}*/

		if (print_mode) {
			graph.streamGraphClustering(cout, bestClustering);
			cout << "Error: " << error << " (" << pError << "/" << nError << ")" << endl;
		}

		if (lastError != -1) {
			if (error > lastError) achievedMinimum = true;
			bool doPrint = false;

			// if (achievedMinimum && hasAllNegativeBonds(graph, bestClustering)) { cout << "Negative bonds on non-optimal" << endl; doPrint = true; }

			// Overall conjecture checking
			if (error < lastError && achievedMinimum) { cout << "Hypothesis bad." << endl; doPrint = true; }

			// First derivatives
			//if (pError < lastPositiveError) { cout << "Hypothesis bad for first derivatives (p)." << endl; doPrint = true; }
			//if (nError > lastNegativeError) { cout << "Hypothesis bad for first derivatives (n)." << endl; doPrint = true; }

			// Reporting
			if (doPrint) graph.streamGraph(cout);
			if (doPrint) graph.streamGraphClustering(cout, bestClustering);
			if (doPrint) cout << "Error: " << error << " (" << pError << "/" << nError << ")" << endl;
		}
		lastError = error;
		//lastPositiveError = pError;
		//lastNegativeError = nError;
	}
}

/*
	Creates a random clustering for the given graph vertices
	V - input graph vertices/nodes
	k - number of clusterings to produce
	generator - random number generator
*/
//couldn't find a random number generator abstract class
vector<vector<int>> getRandomClustering(vector<int> V, unsigned int k, mt19937_64& generator) {
	assert(k <= V.size());

	vector<vector<int>> clustering;
	vector<int> V_left;

	//Assign a random vertex to each cluster to guarantee each has >=1 vertex
	for (size_t i = 0; i < k; i++) {
		clustering.push_back(vector<int>());
		int indexV = generator() % V_left.size();
		clustering[i].push_back(V_left[indexV]);
		V_left[indexV] = V_left[V_left.size() - 1];
		V_left.erase(V_left.end() - 1);
	}

	//Assign the rest of the vertices to random clusters
	for (size_t i = 0; i < V_left.size(); i++) {
		int indexC = generator() % k;
		clustering[indexC].push_back(V_left[i]);
	}

	return clustering;
}

bool isRefinement(const vector<vector<int>>& bigClustering, const vector<vector<int>>& smallClustering) {
	map<int, int> nodeToClusterBig;
	map<int, int> nodeToClusterSmall;

	for (size_t i = 0; i < bigClustering.size(); i++) {
		for (int n : bigClustering[i]) {
			assert(nodeToClusterBig.count(n) == 0);
			nodeToClusterBig[n] = i;
		}
	}
	for (size_t i = 0; i < smallClustering.size(); i++) {
		for (int n : smallClustering[i]) {
			assert(nodeToClusterSmall.count(n) == 0);
			nodeToClusterSmall[n] = i;
		}
	}

	map<int, int> refinement;
	assert(nodeToClusterBig.size() == nodeToClusterSmall.size());
	for (pair<int, int> p : nodeToClusterBig) {
		int n = p.first;
		int bigC = nodeToClusterBig[n];
		int smallC = nodeToClusterSmall[n];

		if (refinement.count(bigC) == 1) {
			if (refinement[bigC] != smallC) return false;
		}
		else {
			refinement[bigC] = smallC;
		}
	}

	return true;
}

bool findRefinement(const vector<vector<vector<vector<int>>>>& clusterings, const vector<int>& errors) {
	int size = errors.size();
	int lastError = -1;

	/*for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			for(int k = 0; k < size; k++){
				if(i == j || j == k || k == i) continue;
				int max_i, mid_i, min_i;

				vector<int> v = {i, j, k};
				sort(v.begin(), v.end());
				min_i = v[0]; mid_i = v[1]; max_i = v[2];

				const vector<vector<vector<int>>>& bigClusterings = clusterings[max_i];
				const vector<vector<vector<int>>>& midClusterings = clusterings[mid_i];
				const vector<vector<vector<int>>>& smallClusterings = clusterings[min_i];

				bool found = false;
				for(const vector<vector<int>>& bigClustering : bigClusterings){
					for(const vector<vector<int>>& midClustering : midClusterings){
						for(const vector<vector<int>>& smallClustering : smallClusterings){
							if(isRefinement(bigClustering, midClustering) &&
								isRefinement(bigClustering, smallClustering) &&
								isRefinement(midClustering, smallClustering)){
								found = true; break;
							}
						}
						if(found) break;
					}
					if(found) break;
				}

				if(!found) return false;
			}
		}
	}*/


	for (int i = 0; i < size; i++) {
		int error = errors[i];

		if (lastError != -1) {
			assert(i > 0);
			const vector<vector<vector<int>>>& bigClusterings = clusterings[i];
			const vector<vector<vector<int>>>& smallClusterings = clusterings[i - 1];

			bool found = false;
			for (const vector<vector<int>>& bigClustering : bigClusterings) {
				for (const vector<vector<int>>& smallClustering : smallClusterings) {
					if (isRefinement(bigClustering, smallClustering)) {
						found = true; break;
					}
				}
				if (found) break;
			}

			if (!found) {
				for (const vector<vector<int>> bigClustering : bigClusterings) {
					for (const vector<int> bcluster : bigClustering) {
						cout << "[";
						for (int n : bcluster) {
							cout << n << ",";
						}
						cout << "] ";
					}
					cout << endl;
				}
				cout << endl << endl;
				for (const vector<vector<int>> smallClustering : smallClusterings) {
					for (const vector<int> scluster : smallClustering) {
						cout << "[";
						for (int n : scluster) {
							cout << n << ",";
						}
						cout << "] ";
					}
					cout << endl;
				}
				return false;
			}
		}

		lastError = error;
	}

	return true;
}

struct parseYoshikoOutputs_result {
	vector<int> numNodes;
	vector<double> refTimes;
	vector<double> refCosts;
	vector<Graph*> graphs;
};
//solved_dir - name of the directory to walk through for solutions
//stopCounter - how many graphs to process, -1 = unlimited
struct parseYoshikoOutputs_result parseYoshikoOutputs(string solved_dir, int stop_counter = -1, bool verbose = false) {
	vector<int> numNodes;
	vector<double> refTimes;
	vector<double> refCosts;
	vector<Graph*> graphs;
	for (auto& p : experimental::filesystem::recursive_directory_iterator(solved_dir)) {
		if (experimental::filesystem::is_directory(p.status())) continue;
		if (p.path().extension() == ".out") continue;

		experimental::filesystem::path p_out(p.path());
		p_out.replace_extension(".out");
		ifstream f_out(p_out);
		string outstring;
		f_out.seekg(0, ios::end);
		outstring.reserve(f_out.tellg());
		f_out.seekg(0, ios::beg);
		outstring.assign(istreambuf_iterator<char>(f_out), istreambuf_iterator<char>());	//Read the entire file

		size_t pos1 = outstring.find("solving (possibly reduced) instances...");
		if (pos1 == string::npos) continue;
		size_t pos2 = outstring.find("real: ", pos1);
		if (pos2 == string::npos) continue;
		size_t pos3 = outstring.find("s", pos2);
		if (pos3 == string::npos) continue;
		string timestr = outstring.substr(pos2 + 6, pos3 - pos2 - 6);
		double time = atof(timestr.c_str());

		pos1 = outstring.find("total cost (data reduction + ILP or heuristic):");
		if (pos1 == string::npos) continue;
		pos2 = outstring.find("\n", pos1);
		if (pos2 == string::npos) continue;
		string coststr = outstring.substr(pos1 + 48, pos2);
		double cost = atof(coststr.c_str());

		experimental::filesystem::path p_graph(p.path());
		p_graph.replace_extension(".cm");
		string str_graph = p_graph.string();
		str_graph.replace(str_graph.find("solved"), 6, "samples/set1");
		ifstream f_graph(str_graph);
		Graph* g = Graph::readJENAStream(f_graph);
		graphs.push_back(g);

		if(verbose)
			cout << "Reference time: " << time << ", graph size: " << g->getNumNodes() << endl;
		numNodes.push_back(g->getNumNodes());
		refTimes.push_back(time);
		refCosts.push_back(cost);

		stop_counter--;
		if (stop_counter == 0) break;
	}

	return {numNodes, refTimes, refCosts, graphs};
}

//Brute force graphs and save them to files
void bruteForceWithSaving(int starting_size = 1, long long starting_index = 0) {
    int size = starting_size;
    long long index = starting_index;
    long long maxIndex = Graph::numberOfPossibleGraphs(size);

    while(true){
        Graph graph(size, index);
        stringstream path;
        path << "saved_graphs/" << size;

        experimental::filesystem::create_directories(path.str());

        path << "/" << index << ".out";
        ofstream f;
        f.open(path.str());

        int lastError = -1;
		int lastPositiveError = -1;
		int lastNegativeError = -1;
        bool achievedMinimum = false;
        vector<vector<vector<vector<int>>>> clusteringsForK;
        vector<int> errorsForK;

		vector<int> pErrors;
		vector<int> nErrors;

        for(int i = 0; i < size; i++){
            vector<vector<int>> bestClustering = bruteForce(graph, i+1);
            clusteringsForK.push_back(bruteForce_bestClusterings);
			graph.streamGraphClustering(f, bestClustering);

			int pError = graph.getPositiveErrorInClustering(bestClustering);
			int nError = graph.getNegativeErrorInClustering(bestClustering);
			int error = pError + nError;

			//Hypothesis testing - now unneeded, leaving it in for archiving reasons
            if(lastError != -1){
                if(error > lastError) achievedMinimum = true;
                if(error < lastError && achievedMinimum) cout << "Hypothesis bad." << size << " " << index << endl;
				if(pError < lastPositiveError) cout << "Hypothesis bad." << size << " " << index << endl;
				if(nError > lastNegativeError) cout << "Hypothesis bad." << size << " " << index << endl;
            }
            lastError = error;
			lastPositiveError = pError;
			lastNegativeError = nError;
            errorsForK.push_back(error);

            f << "Error: " << error << endl;
			pErrors.push_back(pError);
			nErrors.push_back(nError);
        }

        index++;
        if(index >= maxIndex){
            size++;
            maxIndex = Graph::numberOfPossibleGraphs(size);
            index = 0;
        }
    }
}

void bruteForceNonIsomorphic() {
	for (int size = 8; size < 9; size++) {
		stringstream s;
		s << "nonisomorph_graphs/graph" << size << ".g6";
		fstream f(s.str());
		vector<Graph*> graphs = Graph::readG6Stream(f);

		int progress = 0;
		for (Graph* g : graphs) {
			findClustering(*g);
			delete g;

			cout << progress << "/" << graphs.size() << '\r';
			cout.flush();
			progress++;
		}
		cout << size << endl;
	}
}

void gridSearch() {
	string solved_dir = "solved";
	mt19937_64 generator;
	generator.seed(chrono::system_clock::now().time_since_epoch().count());
	Evaluator evaluator;
	vector<string> gen_names; /* = { "GenAlg, 10 pop, 3 kept, 1000 gens, 0.05 mutation",
								"GenAlg, 50 pop, 10 kept, 1000 gens, 0.01 mutation",
								"GenAlg, 50 pop, 3 kept, 10000 gens, 0.10 mutation",
								"GenAlg, 25 pop, 10 kept, 1000 gens, 0.01 mutation" };*/
	vector<function<clustering_t(const Graph&)>> gens;
	vector<int> pops = { 5, 10, 25, 50, 100, 250, 500, 1000 };
	vector<int> kept = { 2, 5, 10, 25, 50 };
	vector<int> geners = { 1, 10, 50, 100, 500 }; //1000, 5000, 10000 };
	vector<double> mutations = { 0.0, 0.0001, 0.001, 0.01, 0.1 };
	GenMemberFactory * factoryNodeList = new GenMemberFactory(GEN_MEMBER_TYPE::NUMBER_LIST);
	GenMemberFactory * factoryCluster = new GenMemberFactory(GEN_MEMBER_TYPE::CLUSTER);
	vector<GenMemberFactory*> factories = { factoryCluster, factoryNodeList };

	for (auto factory : factories) {
		for (int p : pops) {
			for (int k : kept) {
				if (k > p) continue;
				for (int ge : geners) {
					for (double m : mutations) {
						auto gen = [&generator, factory, p, k, ge, m](const Graph& g) {return GenAlgo::geneticAlgorithm(g, generator, factory, p, k, ge, 0.1f, m, false); };
						gens.push_back(gen);
						stringstream ss;
						ss << "GenAlg, " << p << " pop, " << k << " kept, " << ge << " gens, " << m << " mutation";
						gen_names.push_back(ss.str());
					}
				}
			}
		}
	}

	struct parseYoshikoOutputs_result yosh_res = parseYoshikoOutputs(solved_dir, 5);
	
	vector<vector<double>> times; for (size_t i = 0; i < gens.size(); i++) { times.push_back(vector<double>()); }
	vector<vector<double>> relErrors; for (size_t i = 0; i < gens.size(); i++) { relErrors.push_back(vector<double>()); }
	for (size_t i = 0; i < yosh_res.graphs.size(); i++) {
		vector<Graph*> gs = { yosh_res.graphs[i] };

		for (size_t gen_index = 0; gen_index < gens.size(); gen_index++) {
			struct eval_result eval_res = evaluator.evaluateAlgorithm(gens[gen_index], gs);

			double err = eval_res.costs[0];
			double t = eval_res.time;

			double diffErr = err - yosh_res.refCosts[i];

			double relErr = diffErr / yosh_res.refCosts[i];

			times[gen_index].push_back(t);
			relErrors[gen_index].push_back(relErr);

			cout << gen_names[gen_index] << " - time: " << t << ", relErr: " << relErr << endl;
		}
	}

	cout << endl << endl;
	cout << "Graphs and numNodes: [";
	for (int n : yosh_res.numNodes) cout << n << ",";
	cout << "]" << endl;
	cout << "Reference times: [";
	for (double t : yosh_res.refTimes) cout << t << ",";
	cout << "]" << endl;
	for (size_t i = 0; i < gens.size(); i++) {
		cout << gen_names[i] << " times: [";
		for (double t : times[i]) cout << t << ",";
		cout << "]" << endl;

		cout << gen_names[i] << " relErrors: [";
		for (double e : relErrors[i]) cout << e << ",";
		cout << "]" << endl;
	}

	cout << "0 boys -- " << endl;
	for (size_t i = 0; i < gens.size(); i++) {
		bool cont = false;
		for (double e : relErrors[i]) { if (e > 0.5) { cont = true; break; } }
		if (cont) continue;

		cout << gen_names[i] << " times: [";
		for (double t : times[i]) cout << t << ",";
		cout << "]" << endl;
	}
}

int main(int argc, char** argv) {
	string solved_dir = "solved";
	mt19937_64 generator;
	generator.seed(2019);
	Evaluator evaluator;

	GenMemberFactory * factoryNodeList = new GenMemberFactory(GEN_MEMBER_TYPE::NUMBER_LIST);
	GenMemberFactory * factoryCluster = new GenMemberFactory(GEN_MEMBER_TYPE::CLUSTER);
	struct parseYoshikoOutputs_result yoshiko_res = parseYoshikoOutputs(solved_dir, 500);

	vector<function<clustering_t(const Graph&)>> gens;
	vector<string> genNames;
	/*for(float mut : {0.001f, 0.01f, 0.02f, 0.05f, 0.1f}) {
		auto gen = [mut, &generator, factoryNodeList](const Graph& g) {return GenAlgo::geneticAlgorithm(g, generator, factoryNodeList, 250, 2, 50, mut, 0.3f, false); };
		gens.push_back(gen);
		stringstream ss;
		ss << "250 pop, 2 kept, 50 gens, " << mut << " mut, Node list gen";
		genNames.push_back(ss.str());
	}*/
	for (int geners : {100,250,500,1000,2000}) {
		auto gen = [geners, &generator, factoryCluster](const Graph& g) {return GenAlgo::geneticAlgorithm(g, generator, factoryCluster, 5, 2, geners, 0.1f, 0.3f, false); };
		gens.push_back(gen);
		stringstream ss;
		ss << "5 pop, 2 kept, " << geners << " gens, 0.3 mut, Cluster gen";
		genNames.push_back(ss.str());
	}

	struct stats {
		double total_time = 0.0;
		double total_cost = 0.0;
		double total_rel_cost = 0.0;
		double total_rel_cost2 = 0.0;
	};
	vector<struct stats> genStats;
	double total_ref_cost = 0.0;
	double total_ref_time = 0.0;

	for (size_t i = 0; i < yoshiko_res.graphs.size(); i++) {
		total_ref_cost += yoshiko_res.refCosts[i];
		total_ref_time += yoshiko_res.refTimes[i];
	}

	cout << "Total ref cost: " << total_ref_cost << endl;
	cout << "Total ref time: " << total_ref_time << endl;

	for (size_t genIndex = 0; genIndex < gens.size(); genIndex++) {
		struct stats current_stats;

		cout << genNames[genIndex] << ": " << endl;

		for (size_t i = 0; i < yoshiko_res.graphs.size(); i++) {
			vector<Graph*> gs = { yoshiko_res.graphs[i] };
			struct eval_result eval_res = evaluator.evaluateAlgorithm(gens[genIndex], gs);
			assert(eval_res.costs.size() == 1);

			current_stats.total_cost += eval_res.costs[0];

			current_stats.total_rel_cost += (eval_res.costs[0] - yoshiko_res.refCosts[i]) / yoshiko_res.refCosts[i];
			current_stats.total_time += eval_res.time;

			cout << i << "/" << yoshiko_res.graphs.size() << "\r";
			cout.flush();
		}

		current_stats.total_rel_cost2 = (current_stats.total_cost - total_ref_cost) / total_ref_cost;

		cout << "Total cost: " << current_stats.total_cost << endl;
		cout << "Total relative cost: " << current_stats.total_rel_cost << endl;
		cout << "Total total relative cost: " << current_stats.total_rel_cost2 << endl;
		cout << "Total time: " << current_stats.total_time << endl;
	}

	return 0;
}