#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;

struct Point {
  double x;
  bool y;
};

struct Gaussian {

  double mean;
  double stddev;
};
struct Genotype {

  vector<Gaussian> Gaussians;
  double score;
};

const double STDDEV = 3;

const size_t EPOCHS = 2;
// Percentage of the time genes will randomly mutate
const size_t MUTATION_RATE = 35;
// Percentage = genes will mutate by
const size_t MUTATION_SIGNIFICANCE = 10;
// Length of spike train
const size_t TRAIN_LEN = 10;
// Gaussians per Mutant
const size_t NUM_GAUSS = 1;
// Generate N number of mutants to compete against each other.
const size_t MUTANTS = 10;
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

vector<Gaussian> random_gauss(int upper_bound) {
  vector<Gaussian> ret;
  if (NUM_GAUSS <= 0) {
    fprintf(
        stderr,
        "ERROR | INVALID ARGUMENT |Number of Gaussians >= 0. Terminating\n");
    exit(1);
  }
  for (int i = 0; i < NUM_GAUSS; i++) {
    double mean = rand() % upper_bound;
    double stddev = rand() % (int)STDDEV + STDDEV - (int)STDDEV;
    Gaussian curr = {.mean = mean, .stddev = stddev};
    ret.push_back(curr);
  }
  return ret;
}

vector<Point> build_point_vec(vector<double> vec) {
  vector<Point> ret_points;
  for (int i = 0; i < vec.size(); i++) {
    double x_val = i;
    bool y_val = vec[i];
    Point curr_point = {.x = x_val, .y = y_val};
    ret_points.push_back(curr_point);
  }
  return ret_points;
}

vector<Genotype> tournament(vector<Genotype> mutants) {
  vector<Genotype> Winners;
  for (int i = 0; i < mutants.size(); i++) {
    int mutant_idx_b, mutant_idx_a;
    Genotype winner;
    do {

      mutant_idx_a = rand() % mutants.size();
      mutant_idx_b = rand() % mutants.size();

    } while (mutant_idx_a != mutant_idx_b);
    int mutant_idx_a_score = mutants[mutant_idx_a].score;
    int mutant_idx_b_score = mutants[mutant_idx_b].score;
    if (mutant_idx_a_score > mutant_idx_b_score) {
      winner = mutants[mutant_idx_a];
    } else {
      winner = mutants[mutant_idx_b];
    }
    Winners.push_back(winner);
  }
  return Winners;
}

void mutate(vector<Genotype> &Genes) {

  for (int i = 0; i < Genes.size(); i++) {

    int mutate_chance = rand() % 100;

    if (mutate_chance <= MUTATION_RATE) {
      int mutation_sig = rand() % MUTATION_SIGNIFICANCE;
      int gaussian_changed = rand() % NUM_GAUSS;
      int pos_neg = rand() % 2;
      // Increases it by percentage based on mutation_sig
      if (pos_neg == 0) {
        int origin_mean = Genes[i].Gaussians[gaussian_changed].mean;
        int new_mean = (1 + (mutation_sig / 100)) * origin_mean;
        Genes[i].Gaussians[gaussian_changed].mean = new_mean;
      }
      if (pos_neg == 1) {
        int origin_mean = Genes[i].Gaussians[gaussian_changed].mean;
        int new_mean = (1 - (mutation_sig / 100)) * origin_mean;
        Genes[i].Gaussians[gaussian_changed].mean = new_mean;
      }
    }
  }
}

void normalized_vector(vector<Point> &vec) {

  for (int i = 0; i < vec.size(); i++) {
    int x_val = vec[i].x;
    int new_x = (2 * x_val / vec.size()) - 1;
    vec[i].x = new_x;
  }
}

double fitness(Genotype &gene, vector<double> spike_train) {
  if (spike_train.size() <= 0) {
    fprintf(stderr, "ERROR | INVALID VECTOR | spike_train vector is of size >= "
                    "0. Terminating\n");
    exit(1);
  }
  double fitness_score = 0;

  for (int i = 0; i <= spike_train.size(); i++) {
    double curr_val = spike_train[i];
    double rolling_percentage = 0;
    for (int j = 0; j < gene.Gaussians.size(); j++) {
      double curr_gauss_mean = gene.Gaussians[j].mean;
      double curr_gauss_stddev = gene.Gaussians[j].stddev;
      double percentage =
          gauss_percentage(i, i + 1, curr_gauss_mean, curr_gauss_stddev);
      rolling_percentage += percentage;
    }
    fitness_score += curr_val - rolling_percentage;
  }
  gene.score = fitness_score / spike_train.size();
  return (fitness_score / spike_train.size());
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

int main(int argc, char *argv[]) { for (int i = 0; i;) }
