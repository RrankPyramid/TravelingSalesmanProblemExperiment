// Tabu.cpp

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

std::vector<int> generate_initial_solution(int num_cities) {
    std::vector<int> initial_solution(num_cities);
    for (int i = 0; i < num_cities; ++i) {
        initial_solution[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(initial_solution.begin(), initial_solution.end(), g);
    return initial_solution;
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
    int num_cities;
    std::cin >> num_cities;
    std::vector<City> cities(num_cities);
    for (int i = 0; i < num_cities; ++i) {
        double x, y;
        std::cin >> x >> y;
        cities[i] = (City){x, y};
    }

    int max_iterations = 100;
    int tabu_tenure = 50;

    // Generate a random initial solution
    std::vector<int> initial_solution = generate_initial_solution(num_cities);

    // puts("Finish first step");
    // std::cout << "Initial tour length: " << tour_length(cities, initial_solution) << std::endl;
    // std::cout << "First Cities of initial solution: " << std::endl;
    // for(int i = 0; i < 10; i++)
    //     std::cout << initial_solution[i] << " " <<  cities[initial_solution[i]].x << " " << cities[initial_solution[i]].y << std::endl;

    // Apply Tabu Search to the initial solution
    std::vector<int> best_solution = TabuSearch(initial_solution, cities, num_cities, max_iterations, tabu_tenure);

    // puts("Finish second step");

    double best_length = tour_length(cities, best_solution);

    // puts("Finish third step");

    std::cout << "Best tour length: " << best_length << std::endl;
    return 0;

}

