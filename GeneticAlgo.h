#pragma once
#include "graph_utils.h"

#include <random>
#include <vector>
#include <map>
#include <iostream>

class GenMember
{
protected:
	double score;
public:
	GenMember();
	virtual ~GenMember();

	void setScore(double score);
	double getScore() const;
	double getScore(const Graph& graph);
	virtual double computeScore(const Graph & graph) = 0;
	bool operator< (const GenMember &other);
	bool operator> (const GenMember &other);
	virtual clustering_t convertToClustering(const Graph & g) = 0;
};

class GenMemberNumberList : public GenMember {
private:
	int* nodeToCluster;
public:
	GenMemberNumberList(int * nodeToCluster);
	virtual ~GenMemberNumberList();
	static GenMemberNumberList* createRandomMember(const Graph & graph, std::mt19937_64 & generator);
	static GenMemberNumberList* createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMemberNumberList * m1, GenMemberNumberList * m2, float mutation_chance, float crossover_break_chance);
	virtual double computeScore(const Graph & graph);
	virtual clustering_t convertToClustering(const Graph & g);

	int* getNodeToCluster() const;
};

class GenMemberCluster : public GenMember {
private:
	clustering_t clusters;
public:
	GenMemberCluster(clustering_t clusters);
	virtual ~GenMemberCluster();
	clustering_t getClusters();

	static GenMemberCluster* createRandomMember(const Graph & graph, std::mt19937_64 & generator);
	static GenMemberCluster* createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMemberCluster * m1, GenMemberCluster * m2, float mutation_chance, float crossover_break_chance);
	virtual double computeScore(const Graph & graph);
	virtual clustering_t convertToClustering(const Graph & g);
};

enum GEN_MEMBER_TYPE {
	NUMBER_LIST,
	CLUSTER
};

class GenMemberFactory {
private:
	GEN_MEMBER_TYPE member_type;
public:
	GenMemberFactory(GEN_MEMBER_TYPE member_type);
	GenMember * createRandomMember(const Graph & graph, std::mt19937_64 & generator);
	GenMember* createOffspringMember(const Graph & graph, std::mt19937_64 & generator, GenMember * m1, GenMember * m2, float mutation_chance, float crossover_break_chance);
};

class GenAlgo
{
public:
	static clustering_t geneticAlgorithm(const Graph & graph, std::mt19937_64 & generator, GenMemberFactory * factory, size_t num_population, size_t kept_population, size_t num_generations, float stagnation_multiplier, float mutation_rate, bool verbose_mode);
};
