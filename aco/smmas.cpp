#ifndef DHK_SMMAS
#define DHK_SMMAS

#include "aco.cpp"
#include "../standard.cpp"


const double DEFAULT_MIN_PHEROMONE = 0.1;
const double DEFAULT_MAX_PHEROMONE = 1;


struct Smmas : NodeAco{
    double alpha;
    double beta;

    double min_pheromone;
    double max_pheromone;

    Smmas(const vector<int> &layer_sizes,
          const double &rho = DEFAULT_RHO,
          const double &alpha = DEFAULT_ALPHA,
          const double &beta = DEFAULT_BETA,
          const double &min_pheromone = DEFAULT_MIN_PHEROMONE,
          const double &max_pheromone = DEFAULT_MAX_PHEROMONE) :
          NodeAco(layer_sizes, max_pheromone, rho),
          alpha(alpha),
          beta(beta),
          min_pheromone(min_pheromone),
          max_pheromone(max_pheromone) {
    }

    void update_pheromones(const vector<int>& path) override {
        for (int i = 0; i < depth; ++i)
            for (double &pheromone : pheromones[i])
                pheromone += rho*(min_pheromone - pheromone);
        for (int i = 0; i < depth; ++i)
            pheromones[i][path[i]] += rho*(max_pheromone - min_pheromone);
    }

    pair<double, vector<int>> find_path() override {
        // TODO: Implement general finding path process of SMMAS
    }
};


#endif // DHK_SMMAS
