#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;

// Length of spike train
const int TRAIN_LEN = 100;
// Gaussians per Mutant
const int NUM_GAUSS = 5;
// Generate N number of mutants to compete against each other.
const int MUTANTS = 100;
const vector<double> test_data = {0, 0, 1, 0, 1, 1, 1, 0, 1, 0};

// Calculates the bell curve percent for a given region between pos_1 and pos_2
double gauss_percentage(double pos_1, double pos_2, double mean,
                        double stddev) {

  double z_score_1 = (pos_1 - mean) / (stddev);

  double cdf_1 = 0.5 * (1.0 + erf(z_score_1 / sqrt(2)));

  double z_score_2 = (pos_2 - mean) / (stddev);

  double cdf_2 = 0.5 * (1.0 + erf(z_score_2 / sqrt(2)));

  return (cdf_2 - cdf_1);
}

void print_vec(vector<double> vec) {
  printf("| ");
  for (int i = 0; i < vec.size(); i++) {
    printf("%lf ", vec[i]);
  }
  printf("| ");
  printf("\n");
  return;
}
// Score the mutant against the spike train using Mean Absolute Error Algorithm
double score(vector<double> spike_train, vector<double> mutant) {
  double score = 0.0;

  for (size_t i = 0; i < spike_train.size(); ++i) {
    score += abs(spike_train[i] - mutant[i]);
  }

  return 1.0 - (score / spike_train.size());
}

// Generate a N number Gaussians in random locations, summing their percentages.
vector<double> random_vec(int upper_bound) {
  vector<double> ret(upper_bound, 0);

  for (int i = 0; i < NUM_GAUSS; i++) {
    int mean = rand() % upper_bound;
    // Hardcapped stddev to 1 for now.
    int stddev = 1;
    int index_1 = mean;
    int index_2 = index_1 + 1;
    // Runs down the right side of the Gaussian
    while (index_2 <= upper_bound) {
      double result = abs(gauss_percentage(index_1, index_2, mean, stddev));
      ret[index_1] += result;
      index_1++;
      index_2++;
    }
    index_1 = mean - 1;
    index_2 = index_1 - 1;
    // Runs down the left side of the Gaussian
    while (index_2 >= 0) {
      double result = abs(gauss_percentage(index_1, index_2, mean - 1, stddev));
      ret[index_1] += result;
      index_1--;
      index_2--;
    }
    print_vec(ret);
  }
  return ret;
}

int main(int argc, char *argv[]) {

  double best_score = 0;
  vector<double> best_mutant(TRAIN_LEN, 0);
  for (int i = 0; i < MUTANTS; i++) {
    vector<double> curr = random_vec(TRAIN_LEN);
    double scored = score(test_data, curr);
    if (scored > best_score) {
      best_mutant = curr;
      best_score = scored;
    }
  }
  printf("Spike Train: \n");
  print_vec(test_data);
  printf("Best MUTANT: \n");
  print_vec(best_mutant);
  printf("SCORED: %.2lf%% \n", best_score * 100);
}
