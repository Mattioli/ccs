#include <bits/stdc++.h>
#include <chrono>
#include <random>

#define debug(x) std::cerr << #x << " = " << x << std::endl

namespace MtRandomizer {
    int rand(const int a, const int b) {
        static std::seed_seq seq {
            (uint64_t) std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count(),
            (uint64_t) __builtin_ia32_rdtsc()
        };
        static std::mt19937 mt(seq);
        return std::uniform_int_distribution<int>(a, b)(mt);
    }
};


namespace Benchmark {
    template <typename F>
    double benchmark_event(const F&& f) {
        using target_duration = std::chrono::duration<double, std::ratio<1, 1>>;
        
        auto start = std::chrono::steady_clock::now();
        f();
        auto end = std::chrono::steady_clock::now();
        
        return std::chrono::duration_cast<target_duration>(end - start).count();
    }
}

template <typename T>
class Representation {
public:
    Representation() {}

    Representation(const int size) {
        base.assign(size, 0);
    }

    void add(const T index, const T value) {
        if(index < 0 or index >= (int) base.size()) {
            throw std::out_of_range("Out of range access on a Representation object.");
        }
        base[index] += value;
    }

    T size() const {
        return std::accumulate(base.begin(), base.end(), 0);
    }

    void print(std::ostream &out) const {
        for(const T &element : base) {
            out << element << " ";
        }
        out << std::endl;
    }

private:
    std::vector<T> base;
};

template <typename T>
class CoinSystem {
public:
    CoinSystem() {}

    CoinSystem(const std::initializer_list<T> l) {
        base.assign(l);
    }
    
    void insert(const T coin) {
        base.push_back(coin);
        std::sort(base.begin(), base.end());
        base.erase(std::unique(base.begin(), base.end()), base.end());
    }

    Representation<T> greedy_solution(T change) const {
        if(change < 0) {
            throw std::invalid_argument("The change should be a non negative integer.");
        }

        Representation<T> E((int) base.size());
        for(int i = (int) base.size() - 1; i >= 0; i--) {
            T q = change / base[i];
            E.add(i, q);
            change -= q * base[i];
        }

        return E;
    }

    Representation<T> optimal_solution(const T change) const {
        if(change < 0) {
            throw std::invalid_argument("The change should be a non negative integer.");
        }

        std::vector<T> memo(change + 1, -1), parent(change + 1, -1);
        dynamic_solve(change, memo, parent);
        Representation<T> E = construct_solution(change, parent);

        return E;
    }

    bool is_canonical() const {
        if((int) base.size() < 0) {
            throw std::domain_error("The number of coins in the Coin System should be at least 1.");
        }

        if((int) base.size() <= 2) {
            return true;
        }

        T q = base[2] / base[1];
        T r = base[2] % base[1];
        if(0 < r and r < base[2] - q) {
            return false;
        }

        int m = (int) base.size() - 1;
        for(int i = m - 1; i >= 0; i--) {
            for(int j = i; j >= 0; j--) {
                if(base[i] + base[j] > base.back() and greedy_greater_than(base[i] + base[j], 2)) {
                    return false;
                }
            }
        }

        return true;
    }

    bool is_valid() const {
        return base.size() >= 1 && base[0] == 1;
    }

    size_t size() const {
        return base.size();
    }

    void print(std::ostream &out) const {
        for(const T &element : base) {
            out << element << " ";
        }
        out << std::endl;
    }

    static CoinSystem<T> random_system(const int size, const int min_coin_value, const int max_coin_value) {
        CoinSystem<T> C = CoinSystem<T>({1});

        while((int) C.size() < size) {
            int element = MtRandomizer::rand(min_coin_value, max_coin_value);
            C.insert(element);
        }

        return C;
    }

    static CoinSystem<T> random_canonical_system(const int size, const int min_coin_value, const int max_coin_value) {
        CoinSystem<T> C = CoinSystem<T>::random_system(size, min_coin_value, max_coin_value);

        while(!C.is_canonical()) {
            C = CoinSystem<T>::random_system(size, min_coin_value, max_coin_value);
        }

        return C;
    }

private:
    std::vector<T> base;

    bool greedy_greater_than(const T change, const T size) const {
        return greedy_solution(change).size() > size;
    }
    
    T dynamic_solve(T change, std::vector<T> &memo, std::vector<T> &parent) const {
        if(change == 0) {
            return 0;
        }

        T &answer = memo[change];
        if(~answer) {
            return answer;
        }

        answer = std::numeric_limits<T>::max();
        for(int i = 0; i < (int) base.size(); i++) {
            if(change >= base[i]) {
                T current = dynamic_solve(change - base[i], memo, parent) + 1;
                if(current < answer) {
                    answer = current;
                    parent[change] = i;
                }
            }
        }

        return answer;
    }

    Representation<T> construct_solution(T change, const std::vector<T> &parent) const {
        Representation<T> E((int) base.size());

        while(change > 0) {
            T chosen_index = parent[change];
            E.add(chosen_index, 1);
            T chosen_coin = base[chosen_index];
            T next_change = change - chosen_coin;
            change = next_change;
        }

        return E;
    }
};

void generate_bases() {
    const int N_BASES = 50;
    const int MIN_COIN_VALUE = 1;
    const int MAX_COIN_VALUE = (int) 1e5;
    const std::vector<int> ts {5, 10, 25, 50, 100, 250, 500, 1000};
    for(const int &t : ts) {
        freopen(("bases/base_" + std::to_string(t) + ".txt").c_str(), "a+", stdout); 
        std::cout << t << std::endl;
        for(int i = 1; i <= N_BASES; i++) {
            auto base = CoinSystem<long long>::random_canonical_system(t, MIN_COIN_VALUE, MAX_COIN_VALUE);
            base.print(std::cout);
            std::cerr << "Processed base " << i << " of size " << t << "." << std::endl; 
        }
    }
}

int main() {
    const int N_BASES = 50;
    const int N_QUERIES = 500;
    const int MIN_QUERY_VALUE = 0;
    const int MAX_QUERY_VALUE = 10000;
    
    generate_bases();

    for(const int &t : ts) {
        std::vector<CoinSystem<long long>> bases(N_BASES); 
        std::generate(
            bases.begin(),
            bases.end(), 
            std::bind(CoinSystem<long long>::random_canonical_system, t, MIN_COIN_VALUE, MAX_COIN_VALUE)
        );

        std::vector<int> queries(N_QUERIES);
        std::generate(
            queries.begin(),
            queries.end(),
            std::bind(MtRandomizer::rand, MIN_QUERY_VALUE, MAX_QUERY_VALUE)
        );

        for(const auto &base : bases) {
            double dp_total_time = 0.0;
            double greedy_total_time = 0.0;
            double verification_time = Benchmark::benchmark_event(std::bind(&CoinSystem<long long>::is_canonical, base));

            for(const auto &query : queries) {
                double dp_time = Benchmark::benchmark_event(
                    std::bind(&CoinSystem<long long>::optimal_solution, base, query)
                );

                double greedy_time = Benchmark::benchmark_event(
                    std::bind(&CoinSystem<long long>::greedy_solution, base, query)
                );
                
                dp_total_time += dp_time;
                greedy_total_time += greedy_time;
            }

            // TODO: plot A_time and B_time based on t
            double A_time = dp_total_time;
            double B_time = greedy_total_time + verification_time;
        }
    }

    return 0;
}
