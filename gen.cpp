// gen.cpp
// Generates polymers of increasing degree of polymerization (n) 
// and calculates the L_L and L_G values for each polymer
// 
// We switched to C++ for this task for better performance

#include <fstream>
#include <string>
#include <random>
#include <iostream>
#include <getopt.h>

class Args {
private:
    void get_mode(int argc, char * argv[]) {
        opterr = false;

        int choice;
        int index = 0;
        option long_options[] = {
            {"g_prob", required_argument, nullptr, 'g'},
            {"fixed", optional_argument, nullptr, 'f'},
            {"dimers", optional_argument, nullptr, 'd'},
            {nullptr, 0, nullptr, 0}
        };  // long_options[]

        while ((choice = getopt_long(argc, argv, "g:f::d::", long_options, &index)) != -1) {
            switch (choice) {
                case 'g':
                    _g_prob = std::stod(optarg);
                    break;
                case 'f':
                    _fixed = true;
                    if (optarg) {
                        if(optarg[0] == '0' || std::string(optarg) == "false") _fixed = false;
                        else if(optarg[0] == '1' || std::string(optarg) == "true") _fixed = true;
                        else {
                            std::cerr << "Error: invalid option\n";
                            exit(1);
                        }
                    }
                    break;
                case 'd':
                    _dimers = true;
                    if (optarg) {
                        if(optarg[0] == '0' || std::string(optarg) == "false") _dimers = false;
                        else if(optarg[0] == '1' || std::string(optarg) == "true") _dimers = true;
                        else {
                            std::cerr << "Error: invalid option\n";
                            exit(1);
                        }
                    }
                    break;
                case 'h':
                    exit(0);
                default:
                    std::cerr << "Error: invalid option\n";
                    exit(1);
            }  // switch
        }  // while
    }  // getMode()    

    double _g_prob;
    bool _fixed;
    bool _dimers;

public:
    Args(int argc, char * argv[]) {
        _g_prob = 0.25;
        _fixed = false;
        _dimers = false;
        get_mode(argc, argv);
    }  // Args()

    double g_prob() const {
        return _g_prob;
    }  // g_prob()

    bool fixed() const {
        return _fixed;
    }  // fixed()

    bool dimers() const {
        return _dimers;
    }  // dimers()
}; // Args


static std::default_random_engine rng;

// Randomly generate polymer of length N from L and G monomers
// Input: n (int) - length of polymer in monomers (degree of polymerization)
//        g_prob (double) - probability of G monomer occuring at each position
//        fixed (bool) - generate with fixed number of G monomers (fixed vs unfixed method described in paper)
//        dimers (bool) - generate with dimers (true - ring-opening, false - polycondensation)
// Sample runs: (48, 0.25, true, false)  -> LLGLLLGLLLLLGLGLLLLLLLLLLGLLLLLGLGGGGLLGLLLLGLLL
//              (48, 0.25, true, true)   -> LLLLGGLLLLLLLLLLGGLLGGGGLLLLLLLLLLGGLLLLLLLLLLGG
//              (48, 0.25, false, false) -> LLLGGLGLLGLLGLLLLGLLLLLLLLLLLLLGLLGLLLGLLGGGGLLL
std::string gen(int n, 
                double g_prob, 
                bool fixed, 
                bool dimers) {
    std::string polymer;

    if (dimers) n /= 2;
    
    polymer.resize(n, 'L');
    
    if(fixed) {
        std::vector<int> dist(n);
        iota( dist.begin(), dist.end(), 1 );
        std::shuffle(dist.begin(), dist.end(), rng);

        for(int i = 0; i < n * g_prob; ++i) {
            polymer[dist[i]] = 'G';
        } // for
    } else {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for(int i = 0; i < n; ++i) {
            if(dist(rng) < g_prob) {
                polymer[i] = 'G';
            }
        } // for
    } // if...else

    if(!dimers) return polymer;

    std::string final_polymer;
    final_polymer.reserve(n * 2);

    for(int i = 0; i < n; ++i) {
        final_polymer += polymer[i];
        final_polymer += polymer[i];
    } // for

    return final_polymer;
} // gen()

struct Stats {
    int GGs;
    int LLs;
    int GLs;
    int LGs;
}; // Stats

// Calculate GG, LL, GL, and LG counts for a given polymer
// Input: polymer (string) - polymer formed by G and L monomers
Stats calc_stats(const std::string& polymer) {
    Stats stats = {0, 0, 0, 0};
    for(int i = 0; i < polymer.size() - 1; ++i) {
        if(polymer[i] == 'G' && polymer[i + 1] == 'G') {
            stats.GGs++;
        } else if(polymer[i] == 'L' && polymer[i + 1] == 'L') {
            stats.LLs++;
        } else if(polymer[i] == 'G' && polymer[i + 1] == 'L') {
            stats.GLs++;
        } else if(polymer[i] == 'L' && polymer[i + 1] == 'G') {
            stats.LGs++;
        } // if...else
    } // for
    return stats;
} // calc_stats()

double mean(const std::vector<double>& data) {
    double sum = 0;
    for(int i = 0; i < data.size(); ++i) {
        sum += data[i];
    } // for
    return sum / data.size();
} // mean()

double stdev(const std::vector<double>& data, double m) {
    double sum = 0;
    for(int i = 0; i < data.size(); ++i) {
        sum += (data[i] - m) * (data[i] - m);
    } // for
    return sqrt(sum / (data.size()));
} // stdev()

double sem(const std::vector<double>& data, double m) {
    double s = stdev(data, m);
    return s / sqrt(data.size() - 1);
} // sem()

// Calculate L_L or L_G values for a given polymer
// Input: top (vector<int>) - vector of counts of LL or LG
//        bot (vector<int>) - vector of counts of GL or GG
std::vector<double> calc_L_L_or_L_G(std::vector<int>& top, 
                                std::vector<int>& bot) {
    std::vector<double> L_L_or_L_Gs;
    for(int i = 0; i < top.size(); ++i) {
        if(bot[i] == 0) bot[i] = 1;
        L_L_or_L_Gs.push_back((double)top[i] / (double)bot[i] + 1);
    } // for
    return L_L_or_L_Gs;
} // calc_L_L_or_L_G()

int main(int argc, char *argv[]) {
    rng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::ios_base::sync_with_stdio(false);

    Args args(argc, argv);
    int N = 10000;

    std::vector<double> L_L_means;
    std::vector<double> L_L_sems;
    std::vector<double> L_G_means;
    std::vector<double> L_G_sems;

    for(int n = 40; n <= 3000; n += 8) {
        std::vector<int> LL_stats;
        std::vector<int> LG_stats;
        std::vector<int> GG_stats;
        std::vector<int> GL_stats;

        for(int i = 0; i < N; ++i) {
            std::string polymer = gen(n, args.g_prob(), args.fixed(), args.dimers());
            Stats stats = calc_stats(polymer);

            LL_stats.push_back(stats.LLs);
            LG_stats.push_back(stats.LGs);
            GG_stats.push_back(stats.GGs);
            GL_stats.push_back(stats.GLs);
        } // for

        std::vector<double> L_Ls = calc_L_L_or_L_G(LL_stats, LG_stats);
        std::vector<double> L_Gs = calc_L_L_or_L_G(GG_stats, GL_stats);
        
        double L_L_mean = mean(L_Ls);
        L_L_means.push_back(L_L_mean);
        L_L_sems.push_back(sem(L_Ls, L_L_mean));

        double L_G_mean = mean(L_Gs);
        L_G_means.push_back(L_G_mean);
        L_G_sems.push_back(sem(L_Gs, L_G_mean));
    } // for

    std::ofstream file;
    
    std::string append = "";
    if(args.fixed()) append += "_f";
    if(args.dimers()) append += "_d";

    std::cout << L_L_means.size() << std::endl;
    file.open("data/L_L_means" + append + ".txt");
    for(int i = 0; i < L_L_means.size(); ++i) {
        file << L_L_means[i] << "\n";
    } // for
    file.close();

    file.open("data/L_L_sems" + append + ".txt");
    for(int i = 0; i < L_L_sems.size(); ++i) {
        file << L_L_sems[i] << "\n";
    } // for
    file.close();

    file.open("data/L_G_means" + append + ".txt");
    for(int i = 0; i < L_G_means.size(); ++i) {
        file << L_G_means[i] << "\n";
    } // for
    file.close();

    file.open("data/L_G_sems" + append + ".txt");
    for(int i = 0; i < L_G_sems.size(); ++i) {
        file << L_G_sems[i] << "\n";
    } // for
    file.close();
} // main()