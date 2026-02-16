#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;

struct Point {
  double x;
  double y;
};

struct Gaussian {

  double mean;
  double stddev;
  double amp;
};
struct Genotype {

  vector<Gaussian> Gaussians;
  double score;
};

int DENSITY;

double AMP_MAX = 2.5;
double STDDEV = 1;

size_t CROSSBREED_RATE = 75;

size_t EPOCHS = 1000;
// Percentage of the time genes will randomly mutate
size_t MUTATION_RATE = 50;
// Percentage = genes will mutate by
size_t MUTATION_SIGNIFICANCE = 25;
// Length of spike train
size_t TRAIN_LEN = 100;
// Gaussians per Mutant
size_t NUM_GAUSS = 50;
// Generate N number of mutants to compete against each other.
size_t MUTANTS = 1000;

double pdf_calculation(double pos, double mean, double stddev) {

  if (stddev <= 0) {
    return 0;
  }
  double square_stddev = stddev * stddev;
  double coefficient = 1.0 / (stddev * sqrt(2.0 * M_PI));
  double exponent = -((pos - mean) * (pos - mean)) / (2.0 * square_stddev);

  return coefficient * exp(exponent);
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

vector<Genotype> crossbreed(vector<Genotype> generation) {
  vector<Genotype> gen_two;
  for (int i = 0; i + 1 < generation.size(); i += 2) {
    int crossbreed_chance = rand() % 100;

    Genotype parent_a = generation[i];
    Genotype parent_b = generation[i + 1];
    if (crossbreed_chance <= CROSSBREED_RATE) {

      int gauss_to_swap_a = rand() % NUM_GAUSS;
      int gauss_to_swap_b = rand() % NUM_GAUSS;

      Genotype child_a = parent_a;
      child_a.Gaussians[gauss_to_swap_a] = parent_b.Gaussians[gauss_to_swap_b];
      Genotype child_b = parent_b;
      child_b.Gaussians[gauss_to_swap_b] = parent_a.Gaussians[gauss_to_swap_a];
      gen_two.push_back(child_a);
      gen_two.push_back(child_b);
    } else {
      gen_two.push_back(parent_a);
      gen_two.push_back(parent_b);
    }
  }
  return gen_two;
}

vector<Gaussian> random_gauss() {
  vector<Gaussian> ret;
  for (int i = 0; i < NUM_GAUSS; i++) {
    double mean = ((double)rand() / RAND_MAX) * 2.0;
    double stddev = ((double)rand() / RAND_MAX) * STDDEV;
    double amp = ((double)rand() / RAND_MAX) * AMP_MAX;
    ret.push_back({mean, stddev, amp});
  }
  return ret;
}

vector<Point> build_point_vec(vector<double> vec) {
  vector<Point> ret_points;
  for (int i = 0; i < vec.size(); i++) {
    double x_val = i;
    double y_val = vec[i];
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

    } while (mutant_idx_a == mutant_idx_b);
    double mutant_idx_a_score = mutants[mutant_idx_a].score;
    double mutant_idx_b_score = mutants[mutant_idx_b].score;
    if (mutant_idx_a_score < mutant_idx_b_score) {
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
      double delta = mutation_sig / 100.0;
      if (pos_neg == 1) {
        delta = -delta;
      }

      Genes[i].Gaussians[gaussian_changed].mean *= (1.0 + delta);
      Genes[i].Gaussians[gaussian_changed].stddev *= (1.0 + delta);
      Genes[i].Gaussians[gaussian_changed].amp *= (1.0 + delta);
    }
  }
}
void normalized_vector(vector<Point> &vec) {

  for (int i = 0; i < vec.size(); i++) {
    int x_val = vec[i].x;
    double new_x = (2.0 * x_val / vec.size()) - 1;
    vec[i].x = new_x;
  }
}

double fitness(Genotype &gene, const vector<Point> &spike_train) {
  double total_error = 0.0;

  for (size_t i = 0; i < spike_train.size(); i++) {
    double rolling = 0.0;

    for (size_t j = 0; j < gene.Gaussians.size(); j++) {
      rolling += gene.Gaussians[j].amp *
                 pdf_calculation(spike_train[i].x, gene.Gaussians[j].mean,
                                 gene.Gaussians[j].stddev);
    }
    // Normalization of scores
    rolling = 1.0 / (1.0 + exp(-rolling));

    double error = spike_train[i].y - rolling;
    total_error += error * error;
  }

  double score = total_error / spike_train.size();
  gene.score = score;
  return score;
}
int main(int argc, char *argv[]) {

  if (argc != 9) {
    fprintf(
        stderr,
        "Usage: ./Executable <Spike_Train_Length> <Density> <Cross_Breed_Rate> "
        "<Epochs> <Mutation_Rate> <Mutation_Significance> "
        "<Number_of_Gaussians> <Number_of_Mutants>");
  }
  TRAIN_LEN = stoi(argv[1]);
  DENSITY = stoi(argv[2]);

  CROSSBREED_RATE = stoi(argv[3]);

  EPOCHS = stoi(argv[4]);
  MUTATION_RATE = stoi(argv[5]);
  MUTATION_SIGNIFICANCE = stoi(argv[6]);
  NUM_GAUSS = stoi(argv[7]);
  MUTANTS = stoi(argv[8]);

  vector<double> test_data;
  for (int i = 0; i < TRAIN_LEN; i++) {
    int random = rand() % 100;
    if (random <= DENSITY) {
      test_data.push_back(1);
    } else {
      test_data.push_back(0);
    }
  }

  // Spike Train Normalization;
  vector<Point> test_data_p;
  for (int i = 0; i < test_data.size(); i++) {
    test_data_p.push_back({(double)i, test_data[i]});
  }
  normalized_vector(test_data_p);

  vector<Genotype> generation;
  for (int i = 0; i < MUTANTS; i++) {
    vector<Gaussian> gausses = random_gauss();
    Genotype curr = {.Gaussians = gausses, .score = -1};
    fitness(curr, test_data_p);
    generation.push_back(curr);
  }

  for (int i = 0; i < EPOCHS; i++) {
    vector<Genotype> winners = tournament(generation);
    vector<Genotype> crossbred_generation = crossbreed(winners);
    mutate(crossbred_generation);

    double average_fitness = 0;
    double best_fitness = 10000000;
    for (int j = 0; j < winners.size(); j++) {
      double curr_fitness = fitness(crossbred_generation[j], test_data_p);
      average_fitness += curr_fitness;
      if (curr_fitness < best_fitness) {
        best_fitness = curr_fitness;
      }
    }
    average_fitness = average_fitness / winners.size();
    printf("AVERAGE FITNESS: %.2lf\nBEST FITNESS: %.2lf\n", average_fitness,
           best_fitness);

    generation = crossbred_generation;
  }
}
