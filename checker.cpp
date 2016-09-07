#include <algorithm>
#include <cmath>
#include <fstream>
#include <tuple>
#include <vector>

#include "quantile-50-95.hpp"
#include "quantile-50-99.hpp"
#include "quantile-95-95.hpp"

/*
@book{leboudec2010performance,
   title={Performance Evaluation of Computer and Communication Systems},
   author={Le Boudec, Jean-Yves},
   year={2010},
   publisher={EPFL Press, Lausanne, Switzerland}
}
 */
using namespace std;

bool quantile_95_confidence_95(const std::vector<double>& data) {
  /*
  n≤58: no confidence interval possible.

n≥124    ≈⌊0.95n−0.427√n⌋   ≈⌈0.95n+1+0.427√n⌉  0.950
  */
  size_t n = data.size();
  size_t lower_idx;
  size_t upper_idx;

  if (n <= 58)
    return false;
  if (n >= 124) {
    lower_idx = trunc(0.95 * n - 0.427 * sqrt(n));
    upper_idx = trunc(0.95 * n + 1 + 0.427 * sqrt(n));
  }
  const auto it =
      find_if(table_95_95.begin(), table_95_95.end(),
              [n](tuple<int, int, int, double> t) { return get<0>(t) == n; });
  lower_idx = get<1>(*it);
  upper_idx = get<2>(*it);

  printf("%f %f %f\n", data[lower_idx], data[trunc(0.95 * n)], data[upper_idx]);

  return true;
}

bool quantile_50_confidence_95(const std::vector<double>& data) {
  /*
  n≤5: no confidence interval possible.

n≥71    ≈⌊0.50n−0.980√n⌋   ≈⌈0.50n+1+0.980√n⌉  0.950
  */
  size_t n = data.size();
  size_t lower_idx;
  size_t upper_idx;

  if (n <= 5)
    return false;
  if (n >= 71) {
    lower_idx = trunc(0.5 * n - 0.980 * sqrt(n));
    upper_idx = trunc(0.5 * n + 1 + 0.980 * sqrt(n));
  }
  const auto it =
      find_if(table_50_95.begin(), table_50_95.end(),
              [n](tuple<int, int, int, double> t) { return get<0>(t) == n; });
  lower_idx = get<1>(*it);
  upper_idx = get<2>(*it);

  printf("%f %f %f\n", data[lower_idx], data[trunc(n / 2)], data[upper_idx]);

  return true;
}

bool quantile_50_confidence_99(const std::vector<double>& data) {
  /*
  n≤7: no confidence interval possible.

n≥73   ≈⌊0.50n−1.288√n⌋   ≈⌈0.50n+1+1.288√n⌉  0.990
  */
  size_t n = data.size();
  size_t lower_idx;
  size_t upper_idx;

  if (n <= 7)
    return false;
  if (n >= 73) {
    lower_idx = trunc(0.5 * n - 1.288 * sqrt(n))-1;
    upper_idx = trunc(0.5 * n + 1 + 1.288 * sqrt(n))-1;
  }
  const auto it =
      find_if(table_50_99.begin(), table_50_99.end(),
              [n](tuple<int, int, int, double> t) { return get<0>(t) == n; });
  lower_idx = get<1>(*it)-1;
  upper_idx = get<2>(*it)-1;

  printf("%f %f %f\n", data[lower_idx], data[trunc(n / 2)], data[upper_idx]);

  return true;
}

int main(int argc, char** argv) {

  vector<double> data;

  if(argc != 2) {
    printf("usage: checker.x <input file>\n");
    exit(1);
  }
  
  ifstream infile(argv[1]);

  double a;
  while (infile >> a) {
    data.push_back(a);
  }

  sort(data.begin(), data.end());

  quantile_50_confidence_95(data);
  quantile_50_confidence_99(data);
  quantile_95_confidence_95(data);

  return 0;
}
