#ifndef ACOMS_FUNCTION
#define ACOMS_FUNCTION

#include "../standard.cpp"


namespace Function {
    double information_content(const int &n,
                               const int &alp_size,
                               const int &w,
                               const vector<vector<int>> &occurrences,
                               const vector<double> &background,
                               const double &pseudo_count = 0.1) {
        double result = 0.0;

        for (int i = 0; i < alp_size; ++i)
            for (int j = 0; j < w; ++j) {
                double frequency = (double) (occurrences[i][j] + pseudo_count) / (n + pseudo_count*alp_size);
                result += frequency*log(frequency/background[i]);
            }

        return result;
    }
}


#endif // ACOMS_FUNCTION