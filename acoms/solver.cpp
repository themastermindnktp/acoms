#ifndef ACOMS_SOLVER
#define ACOMS_SOLVER

#include "acoms.cpp"
#include "problem.cpp"
#include "../standard.cpp"

namespace Solver {
    const int DEFAULT_N_ANTS = 20;
    const int DEFAULT_N_GENERATIONS = 100;

    string output_filename;

    int n_ants = DEFAULT_N_ANTS;
    int n_generations = DEFAULT_N_GENERATIONS;

    Problem *problem;
    Acoms *acoms;


    void init_from_arguments(int argc, char* argv[]) {
        string dataset_filename;
        int w;

        string alphabet = DNA_ALPHABET;

        double rho = DEFAULT_RHO;
        double alpha = DEFAULT_ALPHA;
        double beta = DEFAULT_BETA;

        double min_pheromone = DEFAULT_MIN_PHEROMONE;
        double max_pheromone = DEFAULT_MAX_PHEROMONE;

        vector<double> background = DEFAULT_DNA_BACKGROUND;

        for (int i = 0; i < argc; ++i) {
            string argument = string(argv[i]);
            if (argument == "-i") dataset_filename = string(argv[++i]);
            else
            if (argument == "-o") output_filename = string(argv[++i]);
            else
            if (argument == "-w") w = atoi(argv[++i]);
            else
            if (argument == "-rho") rho = atof(argv[++i]);
            else
            if (argument == "-alpha") alpha = atof(argv[++i]);
            else
            if (argument == "-beta") beta = atof(argv[++i]);
            else
            if (argument == "-ant") n_ants = atoi(argv[++i]);
            else
            if (argument == "-gen") n_generations = atoi(argv[++i]);
            else
            if (argument == "-min-pher") min_pheromone = atof(argv[++i]);
            else
            if (argument == "-max-pher") max_pheromone = atof(argv[++i]);
            else
            if (argument == "--recalc-background") background.clear();
            else
            if (argument == "--dna") alphabet = DNA_ALPHABET;
        }

        problem = new Problem(dataset_filename, w, alphabet, background);
        acoms = new Acoms(*problem, rho, alpha, beta, min_pheromone, max_pheromone);

        cout << "Processing " << dataset_filename << endl;
        cout << "Number of sequence\t" << problem->n << endl;
        cout << "Total characters  \t" << problem->total << endl;
        cout << "Motif length      \t" << w << endl;
        cout << "Rho               \t" << rho << endl;
        cout << "Alpha             \t" << alpha << endl;
        cout << "Beta              \t" << beta << endl;
        cout << "Min Pheromone     \t" << min_pheromone << endl;
        cout << "Max Pheromone     \t" << max_pheromone << endl;
        cout << endl;

        cout << "Number of ants: " << n_ants << endl;
        cout << "Number of generations: " << n_generations << endl;
        cout << endl;

    }

    void print_result(const vector<vector<int>> &offsets) {
        ofstream output_file;
        if (output_filename.empty()) output_file.copyfmt(cout);
        else {
            cout << "Result saved in " << output_filename << endl;
            output_file.open(output_filename);
        }

        for (int i = 0; i < offsets.size(); ++i) {
            output_file << i + 1 << "=";
            {
                for (int j = 1; j < offsets[i].size(); ++j)
                    output_file << offsets[i][j - 1] << ',';
                output_file << offsets[i].back() << endl;
            }
        }
    }

    void solve() {
        cout << "Running ACOMS:" << endl;

        pair<double, vector<int>> result = acoms->find_best_path(n_ants, n_generations);
        cout << "Information content: " << result.first << endl;

        vector<vector<int>> offsets;
        for (int &offset : result.second) offsets.emplace_back(vector<int>{offset});

        print_result(offsets);
        cout << endl;
    }
}

#endif // ACOMS_SOLVER