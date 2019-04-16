#pragma once
#include "graph_utils.h"
#include <vector>
#include <chrono>
#include <functional>

struct eval_result {
	double time;
	std::vector<double> costs;
};

class Evaluator
{
public:
	Evaluator();
	~Evaluator();

	struct eval_result evaluateAlgorithm(const std::function <clustering_t(const Graph&)>& func, std::vector<Graph*> dataset);
};

