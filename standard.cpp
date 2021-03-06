#ifndef DHK_STANDARD
#define DHK_STANDARD

#include <bits/stdc++.h>

using namespace std;


namespace Standard {
    random_device s_random_device;
    mt19937 s_mt19937(s_random_device());

    double real_random(const double &low, const double &high) {
        uniform_real_distribution<double> dist(low, high);
        return dist(s_mt19937);
    }


    int random_choose(const vector<double> &c_prob) {
        // c_prob: cumulative probabilities whose sum is 1.0
        double pivot = real_random(0.0, 1.0);
        return min((int) (lower_bound(c_prob.begin(), c_prob.end(), pivot) - c_prob.begin()), (int) c_prob.size() - 1);
    }
}


#endif // DHK_STANDARD