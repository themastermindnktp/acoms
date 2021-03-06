#ifndef ACOMS_SOLVER
#define ACOMS_SOLVER

#include "acoms.cpp"
#include "problem.cpp"
#include "../standard.cpp"

const int DEFAULT_N_ANTS = 20;
const int DEFAULT_N_GENERATIONS = 100;
const double DEFAULT_OFFSET_THRESHOLD = 0.01;


namespace Solver {
    string output_filename;

    int n_ants = DEFAULT_N_ANTS;
    int n_generations = DEFAULT_N_GENERATIONS;
    double offset_threshold = DEFAULT_OFFSET_THRESHOLD;

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
            if (argument == "-offset") offset_threshold = atof(argv[++i]);
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
                if (!offsets[i].empty()) output_file << offsets[i].back();
                output_file << endl;
            }
        }
    }

    vector<vector<int>> get_offsets(vector<int> &path) {
        vector<vector<int>> occurrences(problem->alp_size, vector<int>(problem->w));

        for (int i = 0; i < problem->n; ++i)
            for (int k = 0; k < problem->w; ++k)
                occurrences[problem->encode(problem->sequences[i][path[i] + k])][k]++;

        vector<vector<int>> offsets(problem->n, vector<int>());

        double threshold = Function::standard_normal_quantile(1.0 - offset_threshold);

        for (int i = 0; i < problem->n; ++i) {
            string &sequence = problem->sequences[i];

            int m = sequence.length() - problem->w + 1;

            vector<double> scores(m);
            double total = 0.0;

            for (int k = 0; k < problem->w; ++k)
                occurrences[problem->encode(sequence[path[i] + k])][k]--;

            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < problem->w; ++k)
                    occurrences[problem->encode(sequence[j + k])][k]++;

                scores[j] = Function::information_content(problem->n,
                                                          problem->alp_size,
                                                          problem->w,
                                                          occurrences,
                                                          problem->background);

                total += scores[j];

                for (int k = 0; k < problem->w; ++k)
                    occurrences[problem->encode(sequence[j + k])][k]--;
            }

            double mean = total / m;

            double variance = 0.0;
            for (double &score : scores) variance += pow(score - mean, 2);
            variance /= m;

            if (variance == 0.0) for (int j = 0; j < m; ++j) offsets[i].emplace_back(j);
            else {
                for (int j = 0; j < m; ++j)
                    if ((scores[j] - mean) / sqrt(variance) > threshold) offsets[i].emplace_back(j);
            }

            for (int k = 0; k < problem->w; ++k)
                occurrences[problem->encode(sequence[path[i] + k])][k]++;
        }

        return offsets;
    }

    void solve() {
        cout << "Running ACOMS:" << endl;

        pair<double, vector<int>> result = acoms->find_best_path(n_ants, n_generations);
        cout << "Information content: " << result.first << endl;

//        vector<vector<int>> offsets = get_offsets(result.second);

        vector<vector<int>> offsets;
        for (int &offset : result.second) offsets.emplace_back(vector<int>{offset});

        print_result(offsets);
        cout << endl;
    }
}

#endif // ACOMS_SOLVER