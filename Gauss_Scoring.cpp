#include <cstdio>
#include <fstream>
using namespace std;

const int MAX_STDDEV = 3;
const int NUM_GAUSS = 2;

const vector<int> test_data = {0, 0, 1, 0, 1, 1, 1, 0, 1, 0};

void print_vec(vector<int> vec) {
  printf("| ");
  for (int i = 0; i < vec.size(); i++) {
    printf("%d ", vec[i]);
  }
  printf("| ");
  printf("\n");
  return;
}

void mutant_expansion(vector<int> &vec) {

  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] > 1) {
      for (int j = 0; j < vec[i]; j++) {
        if (vec[i + j] == 0)
          vec[i + j] = 1;
        if (vec[i - j] == 0)
          vec[i - j] = 1;
      }
      vec[i] = 1;
    }
  }
}

int score(vector<int> spike_train, vector<int> mutant) {
  int score = 0;
  for (int i = 0; i < mutant.size(); i++) {
    if (mutant[i] == spike_train[i]) {
      score++;
    }
  }
  return score;
}

vector<int> random_vec(int upper_bound) {
  vector<int> ret(upper_bound, 0);
  for (int i = 0; i < NUM_GAUSS; i++) {
    int index = rand() % upper_bound;
    int dev = rand() % MAX_STDDEV + 1;
    ret[index] = dev;
  }
  return ret;
}

int main(int argc, char *argv[]) {

  ifstream in_str(argv[1]);
  int mutants = 1000;
  string density_stream = "DEFAULT";
  getline(in_str, density_stream);
  int best_score = 0;
  vector<int> best_mutant;
  printf("%s\n", density_stream.c_str());
  for (int i = 0; i < mutants; i++) {
    vector<int> curr = random_vec(test_data.size());
    printf("Gaussian Data: \n");
    print_vec(curr);
    printf("Expanded Data: \n");
    mutant_expansion(curr);
    print_vec(curr);
    printf("Actual Data: \n");
    print_vec(test_data);
    int scored = score(test_data, curr);
    printf("score: %d\n", scored);

    if (scored > best_score) {
      best_score = scored;
      best_mutant = curr;
    }
  }
  printf("BEST MUTANT:\n");
  print_vec(best_mutant);
  printf("SCORED: %d\n", best_score);
}
