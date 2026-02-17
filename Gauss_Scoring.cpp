#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <time.h>
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

double density = 0.2;

double amplitude_maximum = 2.5;
double standard_deviation = 1;

double crossover_rate = 0.2;
double mutation_rate = 0.05;
double mutation_magnitude = 0.01;

size_t epochs = 1000;
size_t spike_train_length = 100;

size_t num_gaussians = 8;
size_t population_size = 100;

inline double randf() { return rand() / (double)RAND_MAX; }

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
    printf("[");

    for (size_t i = 0; i < vec.size(); i++) {
        printf("%g", vec[i]);

        if (i != vec.size() - 1) {
            printf(", ");
        }
    }

    printf("]\n");
}

vector<Genotype> crossbreed(const vector<Genotype> &generation) {
    vector<Genotype> gen_two(generation.size());

    for (size_t i = 0; i + 1 < generation.size(); i += 2) {
        Genotype parent_a = generation[i];
        Genotype parent_b = generation[i + 1];

        if (randf() <= crossover_rate) {
            size_t swap_index = rand() % num_gaussians;

            Genotype child_a = parent_a;
            child_a.Gaussians[swap_index] = parent_b.Gaussians[swap_index];

            Genotype child_b = parent_b;
            child_b.Gaussians[swap_index] = parent_a.Gaussians[swap_index];

            gen_two[i] = child_a;
            gen_two[i + 1] = child_b;
        } else {
            gen_two[i] = parent_a;
            gen_two[i + 1] = parent_b;
        }
    }

    return gen_two;
}

vector<Gaussian> random_gauss() {
    vector<Gaussian> ret(num_gaussians);

    for (size_t i = 0; i < num_gaussians; i++) {
        double mean = randf() * 2.0 - 1.0;
        double stddev = randf() * standard_deviation;
        double amp = randf() * amplitude_maximum;

        ret[i] = (Gaussian){
            .mean = mean,
            .stddev = stddev,
            .amp = amp,
        };
    }

    return ret;
}

vector<Point> build_point_vec(const vector<double> &vec) {
    vector<Point> ret_points;

    for (size_t i = 0; i < vec.size(); i++) {
        double x_val = i;
        double y_val = vec[i];
        Point curr_point = {.x = x_val, .y = y_val};

        ret_points[i] = curr_point;
    }

    return ret_points;
}

vector<Genotype> tournament(const vector<Genotype> &mutants) {
    vector<Genotype> Winners(mutants.size());

    for (size_t i = 0; i < 2; i++) {
        Winners[i] = mutants[i];
    }

    for (size_t i = 2; i < mutants.size(); i++) {
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

        Winners[i] = winner;
    }

    return Winners;
}

void mutate(vector<Genotype> &Genes) {
    for (size_t i = 0; i < Genes.size(); i++) {
        for (size_t j = 0; j < Genes[i].Gaussians.size(); j++) {
            if (rand() % 100 < mutation_rate) {
                // Mean
                double perturbation = rand() / ((double)RAND_MAX / 2) - 1.0;

                Genes[i].Gaussians[j].mean += perturbation * mutation_magnitude;
            }

            if (rand() % 100 < mutation_rate) {
                // Standard Deviation
                double perturbation = rand() / ((double)RAND_MAX / 2) - 1.0;

                Genes[i].Gaussians[j].stddev +=
                    perturbation * mutation_magnitude;
            }

            if (rand() % 100 < mutation_rate) {
                // Amp
                double perturbation = rand() / ((double)RAND_MAX / 2) - 1.0;

                Genes[i].Gaussians[j].amp += perturbation * mutation_magnitude;
            }
        }
    }
}

void normalized_vector(vector<Point> &vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        int x_val = vec[i].x;
        double new_x = (2.0 * x_val / vec.size()) - 1.0;

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
    size_t seed = time(nullptr);
    int c;

    while (1) {
        auto option_index = 0;
        static struct option long_options[] = {
            {"spike_train_length", required_argument, 0, 'l'},
            {"spike_train_density", required_argument, 0, 'd'},
            {"epochs", required_argument, 0, 'e'},
            {"population_size", required_argument, 0, 'p'},
            {"mutation_rate", required_argument, 0, 'm'},
            {"mutation_magnitude", required_argument, 0, 'g'},
            {"crossover_rate", required_argument, 0, 'c'},
            {"number_of_gaussians", required_argument, 0, 'n'},
            {"seed", required_argument, 0, 's'},
            {0, 0, 0, 0},
        };

        c = getopt_long_only(argc, argv, "l:d:e:p:m:g:c:n:r:?", long_options,
                             &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
        case 0:
            printf("option %s", long_options[option_index].name);
            if (optarg) {
                printf(" with arg%s", optarg);
            }
            printf("\n");
            break;
        case 'l':
            spike_train_length = stoi(optarg);
            break;
        case 'd':
            density = strtod(optarg, nullptr);
            break;
        case 'e':
            epochs = strtoull(optarg, nullptr, 0);
            break;
        case 'p':
            population_size = strtoull(optarg, nullptr, 0);
            break;
        case 'm':
            mutation_rate = strtod(optarg, nullptr);
            break;
        case 'g':
            mutation_magnitude = stoi(optarg, nullptr);
            break;
        case 'c':
            crossover_rate = strtod(optarg, nullptr);
            break;
        case 'n':
            num_gaussians = stoi(optarg);
            break;
        case 's':
            seed = strtoull(optarg, nullptr, 0);
            srand(seed);
            break;
        case '?':
            break;
        default:
            printf("?? getopt returned character code 0%o ??\n", c);

            return 1;
        }
    }

    srand(seed);

    vector<double> test_data(spike_train_length);
    for (size_t i = 0; i < test_data.size(); i++) {
        test_data[i] = randf() <= density;
    }

    // Spike Train Normalization;
    vector<Point> normalized_spike_train(test_data.size());
    for (size_t i = 0; i < test_data.size(); i++) {
        normalized_spike_train[i] = (Point){
            .x = (double)i,
            .y = test_data[i],
        };
    }

    normalized_vector(normalized_spike_train);

    vector<Genotype> generation(population_size);
    for (size_t i = 0; i < generation.size(); i++) {
        generation[i] = (Genotype){
            .Gaussians = random_gauss(),
            .score = 0,
        };

        fitness(generation[i], normalized_spike_train);
    }

    for (size_t i = 0; i < epochs; i++) {
        sort(generation.begin(), generation.end(),
             [](const Genotype &a, const Genotype &b) {
                 return a.score < b.score;
             });

        vector<Genotype> winners = tournament(generation);
        vector<Genotype> crossbred_generation = crossbreed(winners);
        mutate(crossbred_generation);

        double average_fitness = 0;
        double best_fitness = DBL_MAX;
        for (size_t j = 0; j < winners.size(); j++) {
            double curr_fitness =
                fitness(crossbred_generation[j], normalized_spike_train);
            average_fitness += curr_fitness;

            if (curr_fitness < best_fitness) {
                best_fitness = curr_fitness;
            }
        }

        average_fitness = average_fitness / winners.size();

        printf("epoch: %6zu avgerage_fitness: %10g best_fitness: %10g\n", i,
               average_fitness, best_fitness);

        generation = crossbred_generation;
    }
}
