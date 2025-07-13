// main_test_bench.cpp
#include <cstdlib>

// Forward-declare each sim you want to swap in/out:
int bpm_1sth_only_test();
// int another_sim();

int main() {
    // You can switch which sim to run here:
    return bpm_1sth_only_test();
    // return another_sim();
}