#include "Evaluator.h"



Evaluator::Evaluator()
{
}


Evaluator::~Evaluator()
{
}

struct eval_result Evaluator::evaluateAlgorithm(const std::function <clustering_t(const Graph&)>& func, std::vector<Graph*> dataset)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<double> errors;
	for (Graph* g : dataset) {
		clustering_t algo_clustering = func(*g);
		errors.push_back(g->getErrorInClustering(algo_clustering));
	}
	auto finish = std::chrono::high_resolution_clock::now();
	auto time_in_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
	double time_in_secs = time_in_nano / 1000.0f / 1000.0f / 1000.0f;

	return { time_in_secs, errors };
}
