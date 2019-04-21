#include "GeneticAlgo.h"
#include <cassert>
#include <algorithm>
#include <numeric>

GenMember::GenMember()
{
	score = std::numeric_limits<double>::lowest();
}

GenMember::~GenMember()
{
}

void GenMember::setScore(double score)
{
	this->score = score;
}

double GenMember::getScore() const
{
	return score;
}

double GenMember::getScore(const Graph & graph)
{
	if (score == std::numeric_limits<double>::lowest())
		return computeScore(graph);
	else
		return score;
}

bool GenMember::operator<(const GenMember & other)
{
	return this->score < other.getScore();
}

bool GenMember::operator>(const GenMember & other)
{
	return this->score > other.getScore();
}

// Tries to apply a genetic algorithm to find the best clustering - (finds for all k)
clustering_t GenAlgo::geneticAlgorithm(const Graph& graph, std::mt19937_64& generator, GenMemberFactory * factory, size_t num_population = 50, size_t kept_population = 10, size_t num_generations = 100,
	float stagnation_multiplier = 0.1f, float mutation_rate = -1.0f, bool verbose_mode = false) {
	// Initialise random population
	std::vector<GenMember*> population;
	double bestScore = INFINITY;
	int stagnationLength = 0;
	int graphSize = graph.getNumNodes();

	for (size_t i = 0; i < num_population; i++)
		population.push_back(factory->createRandomMember(graph, generator));

	for (size_t generation = 0; generation < num_generations; generation++) {
		// Evaluate population
		for (size_t i = 0; i < num_population; i++)
			population[i]->getScore(graph);

		std::sort(population.begin(), population.end(), [](GenMember * m1, GenMember * m2){ return *m1 < *m2; });
		if (verbose_mode && generation % 10 == 0) std::cout << "Generation: " << generation << " score: " << population[0]->getScore() << ", current mutation rate: " << (mutation_rate == -1.0f ?
			(1.0f / graphSize * (1.0f + stagnation_multiplier * stagnationLength)) : mutation_rate) << std::endl;

		if (population[0]->getScore() < bestScore) {
			bestScore = population[0]->getScore();
			stagnationLength = 0;
		}
		else {
			stagnationLength++;
		}

		for (size_t i = kept_population; i < num_population; i++) {
			delete population[i];
			// Create child from parents - can have the same parent 'twice' - crossover loses meaning then
			population[i] = factory->createOffspringMember(graph, generator, population[generator() % kept_population], population[generator() % kept_population], mutation_rate == -1.0f ?
				(1.0f / graphSize * (1.0f + stagnation_multiplier * stagnationLength)) : mutation_rate, 0.25f / graphSize);
		}
	}

	clustering_t result(population[0]->convertToClustering(graph));

	// Cleanup
	for (size_t i = 0; i < num_population; i++) {
		delete population[i];
	}

	return result;
}

GenMemberNumberList::GenMemberNumberList(int * nodeToCluster) : nodeToCluster(nodeToCluster)
{
	
}

GenMemberNumberList::~GenMemberNumberList()
{
	delete[] nodeToCluster;
}

GenMemberNumberList * GenMemberNumberList::createRandomMember(const Graph & graph, std::mt19937_64 & generator)
{
	size_t size = graph.getNumNodes();
	int* nodeToCluster = new int[size];
	for (size_t i = 0; i < size; i++) nodeToCluster[i] = generator() % size;
	GenMemberNumberList * member = new GenMemberNumberList(nodeToCluster);
	return member;
}

GenMemberNumberList * GenMemberNumberList::createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMemberNumberList * m1, GenMemberNumberList * m2, float mutation_chance, float crossover_break_chance)
{
	size_t size = graph.getNumNodes();
	int* nodeToCluster = new int[size];

	// Crossover
	int currentParentIndex = generator() % 2 + 1;
	for (size_t i = 0; i < size; i++) {
		nodeToCluster[i] = (currentParentIndex == 1 ? m1 : m2)->getNodeToCluster()[i];
		if (generator() / (float)generator.max() < crossover_break_chance) {
			currentParentIndex = currentParentIndex == 1 ? 2 : 1;
		}
	}

	// Mutation
	for (size_t i = 0; i < size; i++) {
		if (generator() / (float)generator.max() < mutation_chance)
			nodeToCluster[i] = generator() % size;
	}

	return new GenMemberNumberList(nodeToCluster);
}

double GenMemberNumberList::computeScore(const Graph & graph)
{
	score = graph.getErrorInClustering(Graph::convertNodeListToClustering(graph.getNumNodes(), nodeToCluster));
	return score;
}

clustering_t GenMemberNumberList::convertToClustering(const Graph & g)
{
	return Graph::convertNodeListToClustering(g.getNumNodes(), nodeToCluster);
}

int * GenMemberNumberList::getNodeToCluster() const
{
	return nodeToCluster;
}

GenMemberFactory::GenMemberFactory(GEN_MEMBER_TYPE member_type) : member_type(member_type)
{
}

GenMember * GenMemberFactory::createRandomMember(const Graph & graph, std::mt19937_64 & generator)
{
	switch (member_type) {
	case NUMBER_LIST:
		return GenMemberNumberList::createRandomMember(graph, generator);
	case CLUSTER:
		return GenMemberCluster::createRandomMember(graph, generator);
	}
	return nullptr;
}

GenMember * GenMemberFactory::createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMember * m1, GenMember * m2, float mutation_chance, float crossover_break_chance)
{
	switch (member_type) {
	case NUMBER_LIST:
		return (GenMember*)GenMemberNumberList::createOffspringMember(graph, generator, (GenMemberNumberList*)m1, (GenMemberNumberList*)m2, mutation_chance, crossover_break_chance);
	case CLUSTER:
		return (GenMember*)GenMemberCluster::createOffspringMember(graph, generator, (GenMemberCluster*)m1, (GenMemberCluster*)m2, mutation_chance, crossover_break_chance);
	}
	return nullptr;
}

GenMemberCluster::GenMemberCluster(clustering_t clusters) : clusters(clusters)
{
}

GenMemberCluster::~GenMemberCluster()
{
}

clustering_t GenMemberCluster::getClusters()
{
	return clusters;
}

GenMemberCluster * GenMemberCluster::createRandomMember(const Graph & graph, std::mt19937_64 & generator)
{
	size_t size = graph.getNumNodes();
	clustering_t clusters;
	for (size_t i = 0; i < size; i++) {
		clusters.push_back(cluster_t());
	}

	for (size_t i = 0; i < size; i++) {
		clusters[generator() % size].push_back(i);
	}

	GenMemberCluster * member = new GenMemberCluster(clusters);
	return member;
}

GenMemberCluster * GenMemberCluster::createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMemberCluster * m1, GenMemberCluster * m2, float mutation_chance, float crossover_break_chance)
{
	size_t size = graph.getNumNodes();
	clustering_t clusters;

	// Crossover
	/*int currentParentIndex = generator() % 2 + 1;
	for (size_t i = 0; i < size; i++) {
		nodeToCluster[i] = (currentParentIndex == 1 ? m1 : m2)->getNodeToCluster()[i];
		if (generator() / (float)generator.max() < crossover_break_chance) {
			currentParentIndex = currentParentIndex == 1 ? 2 : 1;
		}
	}*/
	double score = 0;
	int parent = generator() % 2;
	if (parent == 0) { clusters = m1->getClusters(); score = m1->getScore(graph); }
	if (parent == 1) { clusters = m2->getClusters(); score = m2->getScore(graph); }

	// Mutation - splitting and coalescing clusters
	for (size_t i = 0; i < size; i++) {
		if (generator() / (float)generator.max() < mutation_chance) {
			if (generator() % 2 == 0) {
				//Coalesce
				size_t otherClusterIndex = generator() % size;
				if (otherClusterIndex == i) continue;

				score += graph.getDeltaInClusteringMerge(clusters[i], clusters[otherClusterIndex]);

				for (int ni : clusters[otherClusterIndex]) { clusters[i].push_back(ni); }
				clusters[otherClusterIndex].clear();
			}
			else {
				//Split
				if (clusters[i].size() < 2) continue;
				size_t otherClusterIndex = 0;
				for (; otherClusterIndex < size; otherClusterIndex++) {
					if (clusters[otherClusterIndex].size() == 0) break;
				}
				assert(otherClusterIndex < size);

				size_t splitIndex = generator() % clusters[i].size();
				for (size_t j = 0; j < splitIndex; j++) {
					clusters[otherClusterIndex].push_back(clusters[i][clusters[i].size()-1]);
					clusters[i].pop_back();
				}

				score += graph.getDeltaInClusteringSplit(clusters[i], clusters[otherClusterIndex]);
			}
		}
	}

	GenMemberCluster * offspring = new GenMemberCluster(clusters);
	offspring->setScore(score);

	return offspring;
}

double GenMemberCluster::computeScore(const Graph & graph)
{
	score = graph.getErrorInClustering(clusters);
	return score;
}

clustering_t GenMemberCluster::convertToClustering(const Graph & g)
{
	return clusters;
}
