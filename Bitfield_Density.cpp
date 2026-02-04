#include <cstdio>
#include <fstream>
#include <vector>
using namespace std;

const vector<int> BITFIELD = {1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1,
                              0, 0, 1, 0, 1, 0, 0, 0, 0, 1};
const vector<int> field = {0, 1, 1, 1, 0};
void print_vec(vector<int> vec) {
  printf("| ");
  for (int i = 0; i < vec.size(); i++) {
    printf("%d ", vec[i]);
  }
  printf("| ");
  printf("\n");
  return;
}

int main(int argc, char *argv[]) {

  if (argc != 3) {
    fprintf(stderr, "Usage: ./Executable <std_deviation> <output_file_name>\n");
    return 1;
  }
  string file_name;
  int std_deviation = atoi(argv[1]);
  file_name = argv[2];
  ofstream str(file_name);
  if (!str.is_open()) {
    fprintf(stderr, "File failed to open.\n");
    return 1;
  }

  int push_back = 0;
  vector<int> averages;
  int i = 0;
  for (i = 0; i < field.size(); i++) {

    int curr_density = field[i];
    int j = 0;
    do {
      if (j != 0) {
        curr_density += field[i + j];
        curr_density += field[i - j];
      }
      j++;
    } while (j <= std_deviation);
    averages.push_back(curr_density);
    push_back++;
  }
  printf("std_deviation: %d\n", std_deviation);
  printf("i: %d\n", i);
  printf("push_back: %d\n", push_back);
  print_vec(field);
  print_vec(averages);
  for (int i = 0; i < averages.size(); i++) {
    str << averages[i];
  }
}
