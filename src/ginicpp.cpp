#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>

using namespace Rcpp;

// Helper function to sort individuals within each group based on income
bool sortFunc(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] < b[0]; // Sorting based on the first element (income)
}

// Struct to hold contribution and its labels
struct Contribution {
  double value; // Contribution value
  int k, h, i, j; // Labels

  Contribution(double val, int k_val, int h_val, int i_val, int j_val)
    : value(val), k(k_val), h(h_val), i(i_val), j(j_val) {}
};

// Main function to calculate components
// [[Rcpp::export]]
List giniDecomposition(NumericMatrix data, bool contrib = false) {
  int N = data.nrow(); // Total number of individuals

  // Using a set to find the number of unique group labels
  std::set<int> uniqueGroups;
  for (int i = 0; i < N; ++i) {
    uniqueGroups.insert(data(i, 1));
  }
  int K = uniqueGroups.size(); // Number of groups

  // Organizing data by groups
  std::vector<std::vector<std::vector<double>>> groups(K);
  for (int i = 0; i < N; ++i) {
    int group = data(i, 1) - 1; // Converting group label to zero-based index
    groups[group].push_back({data(i, 0), data(i, 1), data(i, 2)});
  }

  // Sorting each group by income
  for (auto& group : groups) {
    std::sort(group.begin(), group.end(), sortFunc);
  }

  // Initialize variables for storing the sum of contributions
  double g_w = 0.0; // Sum of within-group contributions
  double g_b = 0.0; // Sum of between-group contributions

  // Vectors to store individual contributions for w and b
  std::vector<Contribution> w_contributions;
  std::vector<Contribution> b_contributions;

  // Preallocate memory
  if (contrib) {
    size_t estimated_size = std::pow(N, 2);
    w_contributions.reserve(estimated_size);
    b_contributions.reserve(estimated_size);
  }

  // Conditionally populate contributions
  for (int k = 0; k < K; ++k) {
    for (int h = 0; h < K; ++h) {
      for (size_t i = 0; i < groups[k].size(); ++i) {
        for (size_t j = 0; j < groups[h].size(); ++j) {
          double p_ij_kh = groups[k][i][2] * groups[h][j][2];
          double abs_diff_i_j_k_h = std::abs(groups[k][i][0] - groups[h][j][0]);

          if (abs_diff_i_j_k_h != 0.0) {
            double abs_diff_i_j_k = std::abs(groups[k][i][0] - groups[k][j][0]);
            double abs_diff_j_k_h = std::abs(groups[k][j][0] - groups[h][j][0]);
            double w = p_ij_kh * abs_diff_i_j_k * (abs_diff_i_j_k_h / (abs_diff_i_j_k + abs_diff_j_k_h));
            double b = p_ij_kh * abs_diff_i_j_k_h - w;
            g_w += w;
            g_b += b;
            if (contrib) {
              w_contributions.push_back(Contribution(w, k, h, i, j));
              b_contributions.push_back(Contribution(b, k, h, i, j));
            }
          }
        }
      }
    }
  }

  // Calculate the average income
  double totalIncome = 0.0;
  for (int i = 0; i < N; ++i) {
    totalIncome += data(i, 0); // Assuming income is the first column
  }
  double avgIncome = totalIncome / N;

  // Calculate the factor (2 * N^2 * average income)
  double factor = 2.0 * N * N * avgIncome;

  // Dividing the results by the factor
  g_w /= factor;
  g_b /= factor;

  // Prepare vectors for DataFrame conversion
  DataFrame w_contrib_df, b_contrib_df;
  if (contrib) {
    std::vector<double> w_values;
    std::vector<int> w_k, w_h, w_i, w_j;
    for (const auto& contrib : w_contributions) {
      w_values.push_back(contrib.value / factor);
      w_k.push_back(contrib.k);
      w_h.push_back(contrib.h);
      w_i.push_back(contrib.i);
      w_j.push_back(contrib.j);
    }
    w_contrib_df = DataFrame::create(
      Named("value") = w_values,
      Named("k") = w_k,
      Named("h") = w_h,
      Named("i") = w_i,
      Named("j") = w_j
    );

    std::vector<double> b_values;
    std::vector<int> b_k, b_h, b_i, b_j;
    for (const auto& contrib : b_contributions) {
      b_values.push_back(contrib.value / factor);
      b_k.push_back(contrib.k);
      b_h.push_back(contrib.h);
      b_i.push_back(contrib.i);
      b_j.push_back(contrib.j);
    }
    b_contrib_df = DataFrame::create(
      Named("value") = b_values,
      Named("k") = b_k,
      Named("h") = b_h,
      Named("i") = b_i,
      Named("j") = b_j
    );

    return List::create(
      Named("g_w") = g_w,
      Named("g_b") = g_b,
      Named("w_contributions") = w_contrib_df,
      Named("b_contributions") = b_contrib_df
    );
  } else {
    return List::create(
      Named("g_w") = g_w,
      Named("g_b") = g_b
    );
  }


}
