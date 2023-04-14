#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

struct City {
    double x;
    double y;
};

double euclidean_distance(const City& a, const City& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

double tour_length(const std::vector<City>& cities, const std::vector<int>& tour) {
    double length = 0;
    for (size_t i = 0; i < tour.size(); ++i) {
        const City& current_city = cities[tour[i]];
        const City& next_city = cities[tour[(i + 1) % tour.size()]];
        length += euclidean_distance(current_city, next_city);
    }
    return length;
}

std::vector<int> generate_initial_solution(const std::vector<std::vector<double>>& pheromones, const std::vector<City>& cities, int num_cities, int starting_city, double alpha, double beta, std::mt19937& rng) {
    std::vector<int> tour = {starting_city};
    std::vector<bool> visited(num_cities, false);
    visited[starting_city] = true;

    for (int step = 1; step < num_cities; ++step) {
        int current_city = tour.back();
        double unvisited_sum_probabilities = 0.0;
        std::vector<double> unvisited_probabilities(num_cities, 0.0);

        for (int next_city = 0; next_city < num_cities; ++next_city) {
            if (!visited[next_city]) {
                double pheromone = std::pow(pheromones[current_city][next_city], alpha);
                double heuristic = std::pow(1.0 / euclidean_distance(cities[current_city], cities[next_city]), beta);
                unvisited_probabilities[next_city] = pheromone * heuristic;
                unvisited_sum_probabilities += pheromone * heuristic;
            }
        }

        std::discrete_distribution<> distribution(unvisited_probabilities.begin(), unvisited_probabilities.end());
        int next_city = distribution(rng);
        visited[next_city] = true;
        tour.push_back(next_city);
    }
    return tour;
}

void update_pheromones(std::vector<std::vector<double>>& pheromones, const std::vector<City>& cities, const std::vector<int>& solution, double evaporation_rate) {
    double distance = 0.0;
    for (size_t i = 0; i < solution.size() - 1; ++i) {
        distance += euclidean_distance(cities[solution[i]], cities[solution[i + 1]]);
    }
    distance += euclidean_distance(cities[solution.back()], cities[solution.front()]);

    double delta_pheromone = (1.0 - evaporation_rate) / distance;
    for (size_t i = 0; i < solution.size() - 1; ++i) {
        pheromones[solution[i]][solution[i + 1]] += delta_pheromone;
        pheromones[solution[i + 1]][solution[i]] += delta_pheromone;
    }
    pheromones[solution.back()][solution.front()] += delta_pheromone;
    pheromones[solution.front()][solution.back()] += delta_pheromone;
}

std::vector<int> TabuSearch(std::vector<int>& solution, const std::vector<City>& cities, int num_cities, int max_iterations, int tabu_tenure) {
    std::vector<std::vector<int>> tabu_list(num_cities, std::vector<int>(num_cities, 0));
    std::vector<int> best_solution = solution;
    double best_length = tour_length(cities, best_solution);
    int num_iterations_without_improvement = 0;

    while (num_iterations_without_improvement < max_iterations) {
        ++num_iterations_without_improvement;

        // std::cout << "Best length: " << best_length << std::endl;
        // std::cout << "Number of iterations without improvement: " << num_iterations_without_improvement << std::endl;

        int city1 = -1;
        int city2 = -1;
        double best_move_length = std::numeric_limits<double>::infinity();

        for (int i = 0; i < num_cities - 1; ++i) {
            for (int j = i + 1; j < num_cities; ++j) {
                std::vector<int> new_solution = solution;
                std::swap(new_solution[i], new_solution[j]);
                double new_length = tour_length(cities, new_solution);

                if (tabu_list[solution[i]][solution[j]] == 0 || (tabu_list[solution[i]][solution[j]] > 0 && new_length < best_length)) {
                    if (new_length < best_move_length) {
                        best_move_length = new_length;
                        city1 = i;
                        city2 = j;
                    }
                }
            }
        }


        if (city1 >= 0 && city2 >= 0) {
            std::swap(solution[city1], solution[city2]);
            tabu_list[solution[city1]][solution[city2]] = tabu_tenure;
            tabu_list[solution[city2]][solution[city1]] = tabu_tenure;
        }

        // Decrement tabu tenure for all entries in the tabu list
        for (int i = 0; i < num_cities; ++i) {
            for (int j = 0; j < num_cities; ++j) {
                if (tabu_list[i][j] > 0) {
                    tabu_list[i][j]--;
                }
            }
        }
        
        // puts("????");
        
        double solution_length = tour_length(cities, solution);

        if (solution_length < best_length) {
            best_solution = solution;
            best_length = solution_length;
            num_iterations_without_improvement = 0;
        }
    }
    return best_solution;
}

int main() {
    std::mt19937 rng(time(0));
    int num_cities;
    std::cin >> num_cities;
    std::vector<City> cities(num_cities);
    for (int i = 0; i < num_cities; ++i) {
        double x, y;
        std::cin >> x >> y;
        cities[i] = {x, y};
    }

    // ACO parameters
    double alpha = 1.0;
    double beta = 5.0;
    double evaporation_rate = 0.5;
    double initial_pheromone = 1.0;
    int num_ants = 10;
    int num_aco_iterations = 10;

    // Tabu Search parameters
    int max_tabu_iterations = 100;
    int tabu_tenure = 50;

    // Initialize pheromones
    std::vector<std::vector<double>> pheromones(num_cities, std::vector<double>(num_cities, initial_pheromone));

    std::vector<int> best_solution;
    double best_length = std::numeric_limits<double>::max();

    for (int iter = 0; iter < num_aco_iterations; ++iter) {
        for (int ant = 0; ant < num_ants; ++ant) {
            int starting_city = rand() % num_cities;
            std::vector<int> initial_solution = generate_initial_solution(pheromones, cities, num_cities, starting_city, alpha, beta, rng);
            std::vector<int> refined_solution = TabuSearch(initial_solution, cities, num_cities, max_tabu_iterations, tabu_tenure);
            double refined_solution_length = tour_length(cities, refined_solution);

            if (refined_solution_length < best_length) {
                best_solution = refined_solution;
                best_length = refined_solution_length;
            }

            // Update pheromones using the refined_solution
            update_pheromones(pheromones, cities, refined_solution, evaporation_rate);
        }
    }
    std::cout << "Best tour length: " << best_length << std::endl;
    return 0;
}

