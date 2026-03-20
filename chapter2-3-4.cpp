#include <iostream>
#include <vector>
#include <omp.h>

int main() {
    int x = 0;
    int nthreads = 0;

    std::vector<int> results(omp_get_max_threads(), -1);

#pragma omp parallel private(x)
    {
        int tid = omp_get_thread_num();

        x = 0;
        x += 1;

        results[tid] = x;

#pragma omp single
        {
            nthreads = omp_get_num_threads();
        }
    }

    for (int i = 0; i < nthreads; ++i) {
        std::cout << "thread " << i << ": x = " << results[i] << std::endl;
    }

    return 0;
}