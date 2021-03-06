#ifndef DHK_ACO
#define DHK_ACO

#include "../standard.cpp"


const double DEFAULT_RHO = 0.5;
const double DEFAULT_ALPHA = 1;
const double DEFAULT_BETA = 1;


struct Aco {
    int depth;

    double rho;

    Aco(const int &depth,
        const double &rho = DEFAULT_RHO) :
        depth(depth),
        rho(rho) {
    }

    virtual void reset() {}

    virtual void update_pheromones(const vector<int>& path) {}

    virtual void improve_result(double &score, vector<int> &path) {}

    virtual pair<double, vector<int>> find_path() {}

    virtual pair<double, vector<int>> find_best_path(const int &n_ants) {
        cout << ".";
        pair<double, vector<int>> candidate = find_path();

        double score = candidate.first;
        vector<int> path = candidate.second;

        for (int i = 1; i < n_ants; ++i) {
            cout << ".";
            candidate = find_path();

            if (candidate.first > score) {
                score = candidate.first;
                path = candidate.second;
            }
        }

        improve_result(score, path);
        return make_pair(score, path);
    }

    virtual pair<double, vector<int>> find_best_path(const int &n_ants, const int &n_generations) {
        reset();
        cout << "> Generation 1";
        pair<double, vector<int>> candidate = find_best_path(n_ants);
        cout << endl;

        double score = candidate.first;
        vector<int> path = candidate.second;

        for(int i = 1; i < n_generations; ++i) {
            cout << "> Generation " << i + 1;
            candidate = find_best_path(n_ants);
            cout << endl;

            if (candidate.first > score) {
                score = candidate.first;
                path = candidate.second;
            }

            update_pheromones(candidate.second);
        }

        return make_pair(score, path);
    }
};


struct NodeAco : Aco {
    double default_pheromone;
    vector<vector<double>> pheromones;

    NodeAco(const vector<int> &layer_sizes,
            const double &default_pheromone = 0.0,
            const double &rho = DEFAULT_RHO) :
            Aco(layer_sizes.size(), rho) {
        pheromones.resize(depth);
        for (int i = 0; i < depth; ++i) pheromones[i].resize(layer_sizes[i], default_pheromone);
    }

    void reset() override {
        for (int i = 0; i < depth; ++i)
            for (double &pheromone : pheromones[i])
                pheromone = default_pheromone;
    }
};


#endif // DHK_ACO
