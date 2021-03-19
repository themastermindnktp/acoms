#ifndef ACOMS_FUNCTION
#define ACOMS_FUNCTION


#include "../standard.cpp"

namespace Function {
    double information_content(const int &n, // number of motif instances
                               const int &alp_size,
                               const int &w,
                               const vector<vector<int>> &occurrences,
                               const vector<double> &background,
                               const double &pseudo_count = 1) {
        double result = 0.0;

        for (int i = 0; i < alp_size; ++i)
            for (int j = 0; j < w; ++j) {
                double frequency = (double) (occurrences[i][j] + pseudo_count) / (n + pseudo_count*alp_size);
                result += frequency*log(frequency/background[i]);
            }

        return result;
    }

    // quantile function of standard normal distribution
    double standard_normal_quantile(const double &p) {
        // Beasley-Springer-Moro Algorithm
        static double a[4] = {2.50662823884,
                              -18.61500062529,
                              41.39119773534,
                              -25.44106049637};

        static double b[4] = {-8.47351093090,
                              23.08336743743,
                              -21.06224101826,
                              3.13082909833};

        static double c[9] = {0.3374754822726147,
                              0.9761690190917186,
                              0.1607979714918209,
                              0.0276438810333863,
                              0.0038405729373609,
                              0.0003951896511919,
                              0.0000321767881768,
                              0.0000002888167364,
                              0.0000003960315187};

        double t = p - 0.5;
        if (abs(t) < 0.42) {
            double r = t * t;
            double x = (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) * t;
            double y = (((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1;
            return x / y;
        }
        else {
            double r = t > 0 ? 1 - p : p;
            r = log(-log(r));
            double x = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));
            return (t > 0) ? x : -x;
        }
    }

    const int EXPO = 3;
    const int SIDE = 6;

    vector<double> motif_learning_coefficients;
    vector<double> regularization_center;

    void import_motif_learning(const string &coefficient_file_name,
                               const string &regularization_file_name) {
        ifstream coefficient_file(coefficient_file_name);
        double coefficient;
        while (coefficient_file >> coefficient)
            motif_learning_coefficients.emplace_back(coefficient);

        ifstream regularization_file(regularization_file_name);
        double coordinate;
        while (regularization_file >> coordinate)
            regularization_center.emplace_back(coordinate);
    }

    double motif_learning_score(const int &n,
                                const int &alp_size,
                                const int &w,
                                const vector<vector<int>> &occurrences,
                                const vector<double> &background,
                                const double &pseudo_count = 0) {
        vector<vector<double>> frequencies(alp_size, vector<double>());

        vector<double> parameters {1};

        for (int i = 0; i < alp_size; ++i) {
            for (int j = 0; j < w; ++j)
                frequencies[i].emplace_back((double) (occurrences[i][j] + pseudo_count) / (n + pseudo_count * alp_size));

            for (int j = 0; j < SIDE; ++j) {
                double front_value = frequencies[i][j] / background[i];
                for (int expo = 1; expo <= EXPO; ++expo)
                    parameters.emplace_back(pow(front_value, expo));

                double back_value = frequencies[i][w - j - 1] / background[i];
                for (int expo = 1; expo <= EXPO; ++expo)
                    parameters.emplace_back(pow(back_value, expo));
            }

            for (int j = 0; j < SIDE; ++j) {
                double average_frequency = 0;
                for (int k = j; k <= j + w - SIDE; ++k)
                    average_frequency += frequencies[i][k] / background[i];
                average_frequency /= w - SIDE + 1;
                for (int expo = 1; expo <= EXPO; ++expo)
                    parameters.emplace_back(pow(average_frequency, expo));
            }
        }

        double score = 0;
        for (int i = 0; i < parameters.size(); ++i)
            score += parameters[i]*motif_learning_coefficients[i] - sqrt(pow(parameters[i] - regularization_center[i], 2));

        return score;
    }
}

#endif // ACOMS_FUNCTION