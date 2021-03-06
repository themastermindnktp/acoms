#include "acoms/solver.cpp"

int main(int argc, char* argv[]) {
    Solver::init_from_arguments(argc, argv);
    Solver::solve();
    return 0;
}
