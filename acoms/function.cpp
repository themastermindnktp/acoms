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
}


#endif // ACOMS_FUNCTION