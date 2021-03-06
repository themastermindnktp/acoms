#ifndef ACOMS_ACOMS
#define ACOMS_ACOMS

#include "function.cpp"
#include "problem.cpp"
#include "../standard.cpp"
#include "../aco/smmas.cpp"


struct Acoms : Smmas {
    Problem &problem;

    Acoms(Problem &problem,
          const double &rho = DEFAULT_RHO,
          const double &alpha = DEFAULT_ALPHA,
          const double &beta = DEFAULT_BETA,
          const double &min_pheromone = DEFAULT_MIN_PHEROMONE,
          const double &max_pheromone = DEFAULT_MAX_PHEROMONE) :
          problem(problem),
          Smmas(problem.sequence_lengths(), rho, alpha, beta, min_pheromone, max_pheromone) {
    }

    pair<double, vector<int>> find_path() override {
        vector<vector<int>> occurrences(problem.alp_size, vector<int>(problem.w));

        for (int i = 0; i < problem.alp_size; ++i)
            for (int j = 0; j < problem.w; ++j)
                occurrences[i][j] = 0;

        vector<int> path;

        vector<double> options;

        for (int i = 0; i < problem.n; ++i) {
            string &sequence = problem.sequences[i];
            options.clear();

            double total_value = 0.0;

            for (int j = 0; j <= sequence.length() - problem.w; ++j) { // motif starts from j
                for (int k = 0; k < problem.w; ++k) // position k in the motif is position j+k in the sequence
                    occurrences[problem.encode(sequence[j + k])][k]++;

                double heuristic = Function::information_content(i + 1,
                                                                 problem.alp_size,
                                                                 problem.w,
                                                                 occurrences,
                                                                 problem.background);

                double value = pow(pheromones[i][j], alpha) * pow(heuristic, beta);
                options.emplace_back(value);
                total_value += value;

                for (int k = 0; k < problem.w; ++k)
                    occurrences[problem.encode(sequence[j + k])][k]--;
            }

            options.front() /= total_value;
            for (int j = 1; j < options.size(); ++j)
                options[j] = options[j] / total_value + options[j - 1];

            int chosen_index = Standard::random_choose(options);

            path.emplace_back(chosen_index);

            for (int k = 0; k < problem.w; ++k)
                occurrences[problem.encode(sequence[chosen_index + k])][k]++;

        }

        double score = Function::information_content(problem.n,
                                                     problem.alp_size,
                                                     problem.w,
                                                     occurrences,
                                                     problem.background);

        return make_pair(score, path);
    }
};


#endif //ACOMS_ACOMS