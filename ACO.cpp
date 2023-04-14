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

std::vector<int> ACO_construct_solution(const std::vector<std::vector<double>>& pheromones, const std::vector<City>& cities, int num_cities, int starting_city, double alpha, double beta, std::mt19937& rng) {
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
    // Adjust parameters based on input size
    double alpha = 1.0;
    double beta = 5.0;
    double evaporation_rate = 0.3;
    double initial_pheromone = 1.0;
    int num_ants = 10;
    int num_iterations = 10;

    // Initialize pheromones
    std::vector<std::vector<double>> pheromones(num_cities, std::vector<double>(num_cities, initial_pheromone));

    std::vector<int> best_solution;
    double best_distance = std::numeric_limits<double>::max();

    for (int iter = 0; iter < num_iterations; ++iter) {
        for (int ant = 0; ant < num_ants; ++ant) {
            int starting_city = rng() % num_cities;
            std::vector<int> solution = ACO_construct_solution(pheromones, cities, num_cities, starting_city, alpha, beta, rng);

            double distance = 0.0;
            for (size_t i = 0; i < solution.size() - 1; ++i) {
                distance += euclidean_distance(cities[solution[i]], cities[solution[i + 1]]);
            }
            distance += euclidean_distance(cities[solution.back()], cities[solution.front()]);

            if (distance < best_distance) {
                best_distance = distance;
                best_solution = solution;
            }

            update_pheromones(pheromones, cities, solution, evaporation_rate);
        }

        for (int i = 0; i < num_cities; ++i) {
            for (int j = 0; j < num_cities; ++j) {
                if (i != j) {
                    pheromones[i][j] *= (1.0 - evaporation_rate);
                }
            }
        }
    }

    std::cout << "Best tour length: " << best_distance << std::endl;
    return 0;
}