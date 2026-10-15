// 3-way Circuit Partitioning
// Strategy: Pure SA (< 50000 nets) & FM + global SA (>= 50000 nets)
#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <chrono>
#include <numeric>
#include <random>
#include <algorithm>
#include <array>

using namespace std;

// Strategy 1: Pure SA (For Nets < 50000)
namespace SmallSolver {
    struct Net {
        string name;
        vector<int> cells;
        int dist[3] = {0, 0, 0}; 
    };

    struct Cell {
        string name;
        vector<int> nets;
        int partition = -1;
    };

    double r;
    vector<Net> nets;
    vector<Cell> cells;
    unordered_map<string, int> cell_name_to_id;
    int partition_sizes[3] = {0, 0, 0};
    int total_cut_size = 0;

    int best_cut_size = 1e9;
    vector<int> best_partitions;

    // Check if a net spans across multiple partitions
    inline int get_cut(int n_id) {
        int span = (nets[n_id].dist[0] > 0) + (nets[n_id].dist[1] > 0) + (nets[n_id].dist[2] > 0);
        return span > 1 ? 1 : 0;
    }

    // Evaluate cut size delta for moving a single cell
    int evaluate_single_move(int c, int to_p, bool apply = false) {
        int from_p = cells[c].partition;
        int delta = 0;
        
        for (int n_id : cells[c].nets) {
            int cut_before = get_cut(n_id);
            nets[n_id].dist[from_p]--;
            nets[n_id].dist[to_p]++;
            int cut_after = get_cut(n_id);
            delta += (cut_after - cut_before);
            
            if (!apply) { 
                nets[n_id].dist[from_p]++;
                nets[n_id].dist[to_p]--;
            }
        }
        
        if (apply) {
            cells[c].partition = to_p;
            partition_sizes[from_p]--;
            partition_sizes[to_p]++;
            total_cut_size += delta;
        }
        return delta;
    }

    // Evaluate cut size delta for swapping two cells
    int evaluate_swap(int c1, int c2, bool apply = false) {
        int p1 = cells[c1].partition;
        int p2 = cells[c2].partition;
        int delta = 0;

        vector<int> affected;
        affected.reserve(cells[c1].nets.size() + cells[c2].nets.size());
        affected.insert(affected.end(), cells[c1].nets.begin(), cells[c1].nets.end());
        affected.insert(affected.end(), cells[c2].nets.begin(), cells[c2].nets.end());
        sort(affected.begin(), affected.end());
        affected.erase(unique(affected.begin(), affected.end()), affected.end());

        for (int n_id : affected) delta -= get_cut(n_id);

        for (int n_id : cells[c1].nets) { nets[n_id].dist[p1]--; nets[n_id].dist[p2]++; }
        for (int n_id : cells[c2].nets) { nets[n_id].dist[p2]--; nets[n_id].dist[p1]++; }

        for (int n_id : affected) delta += get_cut(n_id);

        if (!apply) { 
            for (int n_id : cells[c1].nets) { nets[n_id].dist[p1]++; nets[n_id].dist[p2]--; }
            for (int n_id : cells[c2].nets) { nets[n_id].dist[p2]++; nets[n_id].dist[p1]--; }
        } else {
            cells[c1].partition = p2;
            cells[c2].partition = p1;
            total_cut_size += delta;
        }
        return delta;
    }

    int run(int argc, char* argv[]) {
        auto start_time = chrono::steady_clock::now();

        ifstream fin(argv[1]);
        fin >> r; 
        string token;
        while (fin >> token) {
            if (token == "NET") { 
                Net net;
                fin >> net.name;
                while (fin >> token && token != ";") { 
                    int c_id;
                    if (cell_name_to_id.find(token) == cell_name_to_id.end()) {
                        c_id = cells.size();
                        cell_name_to_id[token] = c_id;
                        cells.push_back({token, {}, -1});
                    } else {
                        c_id = cell_name_to_id[token];
                    }
                    net.cells.push_back(c_id);
                }
                sort(net.cells.begin(), net.cells.end());
                net.cells.erase(unique(net.cells.begin(), net.cells.end()), net.cells.end());
                if (net.cells.size() > 1) {
                    nets.push_back(net);
                    int n_id = nets.size() - 1;
                    for (int c_id : net.cells) {
                        cells[c_id].nets.push_back(n_id);
                    }
                }
            }
        }
        fin.close();

        int n = cells.size();
        if (n == 0) return 0;

        double min_limit = (1.0 - r) / 3.0 * n;
        double max_limit = (1.0 + r) / 3.0 * n;
        int lower_bound = ceil(min_limit - 1e-6);
        int upper_bound = floor(max_limit + 1e-6);

        auto is_balanced = [&](int p, int diff) {
            int new_size = partition_sizes[p] + diff;
            return new_size >= lower_bound && new_size <= upper_bound;
        };

        // Initial balanced partitioning
        for (int i = 0; i < n; ++i) {
            int p = i % 3;
            cells[i].partition = p;
            partition_sizes[p]++;
        }

        for (int i = 0; i < (int)nets.size(); ++i) {
            for (int c_id : nets[i].cells) {
                nets[i].dist[cells[c_id].partition]++;
            }
            total_cut_size += get_cut(i);
        }
        
        best_cut_size = total_cut_size;
        best_partitions.resize(n);
        for (int i = 0; i < n; ++i) best_partitions[i] = cells[i].partition;

        mt19937 rng(1337); 
        uniform_int_distribution<int> rand_cell(0, n - 1);
        uniform_int_distribution<int> rand_part(0, 2);
        uniform_real_distribution<double> rand_prob(0.0, 1.0);

        double T = 10.0;             
        double alpha = 0.9995;       
        int moves_per_T = max(1000, 20 * n); 

        int total_iters = 0;
        bool time_limit_reached = false;

        // Simulated Annealing
        while (T > 0.001 && !time_limit_reached) {
            for (int step = 0; step < moves_per_T; ++step) {
                total_iters++;
                
                // Safe exit condition
                if (total_iters % 50000 == 0) {
                    auto curr_time = chrono::steady_clock::now();
                    chrono::duration<double> elapsed = curr_time - start_time;
                    if (elapsed.count() > 3400.0) { 
                        time_limit_reached = true;
                        break;
                    }
                }

                int move_type = uniform_int_distribution<int>(0, 1)(rng);
                
                if (move_type == 0) { 
                    int c = rand_cell(rng);
                    int from_p = cells[c].partition;
                    int to_p = rand_part(rng);
                    
                    if (from_p == to_p) continue;
                    if (!is_balanced(from_p, -1) || !is_balanced(to_p, 1)) continue;

                    int delta = evaluate_single_move(c, to_p, false);

                    if (delta < 0 || rand_prob(rng) < exp(-delta / T)) {
                        evaluate_single_move(c, to_p, true); 
                        if (total_cut_size < best_cut_size) {
                            best_cut_size = total_cut_size;
                            for (int i = 0; i < n; ++i) best_partitions[i] = cells[i].partition;
                        }
                    }
                } else { 
                    int c1 = rand_cell(rng);
                    int c2 = rand_cell(rng);
                    int p1 = cells[c1].partition;
                    int p2 = cells[c2].partition;
                    
                    if (p1 == p2) continue;

                    int delta = evaluate_swap(c1, c2, false);

                    if (delta < 0 || rand_prob(rng) < exp(-delta / T)) {
                        evaluate_swap(c1, c2, true); 
                        if (total_cut_size < best_cut_size) {
                            best_cut_size = total_cut_size;
                            for (int i = 0; i < n; ++i) best_partitions[i] = cells[i].partition;
                        }
                    }
                }
            }
            T *= alpha; 
        }

        ofstream fout(argv[2]);
        fout << "Cutsize = " << best_cut_size << "\n";
        
        vector<vector<string>> ans(3);
        for (int i = 0; i < n; ++i) {
            ans[best_partitions[i]].push_back(cells[i].name);
        }
        
        for (int p = 0; p < 3; ++p) {
            fout << "G" << p + 1 << " " << ans[p].size() << "\n";
            for (size_t i = 0; i < ans[p].size(); ++i) {
                fout << ans[p][i];
                if (i + 1 == ans[p].size()) fout << ";\n";
                else fout << " ";
            }
        }
        
        fout.close();
        return 0;
    }
}

// Strategy 2: 3-way FM Partitioning with ILS & SA (For Nets >= 50000)
namespace LargeSolver {
    struct Net {
        string name;
        vector<int> cells;
    };

    struct Cell {
        string name;
        vector<int> nets;
    };

    struct Move { 
        int cell, from, to, cost_change; 
    };

    int run(int argc, char* argv[]) {
        ifstream fin(argv[1]);
        ofstream fout(argv[2]);

        double r;
        if (!(fin >> r)) return 0;

        unordered_map<string, int> cell_name_to_id;
        vector<Cell> cells;
        vector<Net> nets;

        auto get_cell_id = [&](const string& name) {
            if(cell_name_to_id.find(name) == cell_name_to_id.end()) {
                cell_name_to_id[name] = cells.size();
                cells.push_back({name, {}});
            }
            return cell_name_to_id[name];
        };

        string token;
        while(fin >> token) {
            if(token == "NET") {
                string net_name;
                fin >> net_name;
                nets.push_back({net_name, {}});
                int net_id = (int)nets.size() - 1;
                
                while(fin >> token) {
                    if(token == ";") break;
                    if(!token.empty() && token.back() == ';') {
                        token.pop_back();
                        if(!token.empty()) nets[net_id].cells.push_back(get_cell_id(token));
                        break;
                    }
                    if(!token.empty()) nets[net_id].cells.push_back(get_cell_id(token));
                }
            }
        }

        int n = cells.size();
        if(n == 0) return 0;

        // Remove duplicate cells in nets and discard single-cell nets
        int max_degree = 0;
        int valid_nets = 0;
        for(int i = 0; i < (int)nets.size(); ++i) {
            sort(nets[i].cells.begin(), nets[i].cells.end());
            nets[i].cells.erase(unique(nets[i].cells.begin(), nets[i].cells.end()), nets[i].cells.end());
            
            if (nets[i].cells.size() > 1) { 
                if (valid_nets != i) nets[valid_nets] = std::move(nets[i]);
                for(int c : nets[valid_nets].cells) {
                    cells[c].nets.push_back(valid_nets);
                }
                valid_nets++;
            }
        }
        nets.resize(valid_nets);
        int m = nets.size();

        for(int i = 0; i < n; ++i) {
            max_degree = max(max_degree, (int)cells[i].nets.size());
        }

        int lower_bound = ceil(n * (1.0 - r) / 3.0);
        int upper_bound = floor(n * (1.0 + r) / 3.0);
        if(lower_bound > upper_bound) {
            lower_bound = n / 3;
            upper_bound = (n + 2) / 3;
        }

        vector<int> current_part(n, -1);
        vector<int> best_overall_part(n, 0);
        int best_overall_cut = 1e9;
        
        vector<array<int, 3>> net_counts(m);
        vector<array<int, 3>> cell_gain(n);
        vector<int> part_sizes(3);
        vector<int> net_order(m);
        iota(net_order.begin(), net_order.end(), 0);
        
        vector<bool> locked(n);
        vector<Move> history;
        history.reserve(n); 

        int OFFSET = max_degree * 2;
        int MAX_GAIN = max_degree * 4 + 1;
        vector<int> head(MAX_GAIN);
        vector<int> prev_node(n * 3);
        vector<int> next_node(n * 3);
        vector<uint8_t> in_bucket(n * 3);

        mt19937 rng(1337);
        auto start_time = chrono::steady_clock::now();
        double time_limit = 3400.0; 

        bool force_restart = true;
        vector<int> accepted_part(n, -1);
        int accepted_cost = 1e9;
        
        double T = 100.0;            
        double alpha = 0.95;         
        int non_improve_iters = 0;   
        uniform_real_distribution<double> dist_prob(0.0, 1.0);
        int seed_count = 0;

        while(true) {
            if(chrono::duration<double>(chrono::steady_clock::now() - start_time).count() > time_limit) break;
            seed_count++;

            int current_cost = 0;

            if (force_restart) {
                // Initialize randomly within constraints
                fill(current_part.begin(), current_part.end(), -1);
                fill(part_sizes.begin(), part_sizes.end(), 0);
                shuffle(net_order.begin(), net_order.end(), rng);

                int curr_p = 0;
                for(int nid : net_order) {
                    for(int c : nets[nid].cells) {
                        if(current_part[c] == -1) {
                            if(part_sizes[curr_p] < upper_bound) {
                                current_part[c] = curr_p;
                                part_sizes[curr_p]++;
                            } else {
                                for(int p = 0; p < 3; ++p) {
                                    if(part_sizes[p] < upper_bound) {
                                        current_part[c] = p;
                                        part_sizes[p]++;
                                        curr_p = p;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                for(int i = 0; i < n; ++i) {
                    if(current_part[i] == -1) {
                        for(int p = 0; p < 3; ++p) {
                            if(part_sizes[p] < upper_bound) {
                                current_part[i] = p;
                                part_sizes[p]++;
                                break;
                            }
                        }
                    }
                }

                fill(net_counts.begin(), net_counts.end(), array<int, 3>{0, 0, 0});
                for(int i = 0; i < m; ++i) {
                    for(int c : nets[i].cells) net_counts[i][current_part[c]]++;
                    int span = (net_counts[i][0] > 0) + (net_counts[i][1] > 0) + (net_counts[i][2] > 0);
                    if(span > 1) current_cost++; 
                }
                
                accepted_part = current_part;
                accepted_cost = current_cost;
                force_restart = false;
            } else {
                // Perturb previous state to escape local minima
                current_part = accepted_part;
                fill(part_sizes.begin(), part_sizes.end(), 0);
                fill(net_counts.begin(), net_counts.end(), array<int, 3>{0, 0, 0});
                
                for(int i = 0; i < n; ++i) part_sizes[current_part[i]]++;
                for(int i = 0; i < m; ++i) {
                    for(int c : nets[i].cells) net_counts[i][current_part[c]]++;
                }

                int perturb_cnt = max(10, n / 100);
                uniform_int_distribution<int> cell_dist(0, n - 1);
                uniform_int_distribution<int> part_dist(0, 2);

                for(int i = 0; i < perturb_cnt; ++i) {
                    int c = cell_dist(rng);
                    int old_p = current_part[c];
                    int new_p = part_dist(rng);
                    if(old_p == new_p) new_p = (new_p + 1) % 3;

                    if(part_sizes[old_p] - 1 >= lower_bound && part_sizes[new_p] + 1 <= upper_bound) {
                        current_part[c] = new_p;
                        part_sizes[old_p]--;
                        part_sizes[new_p]++;
                        for(int net_id : cells[c].nets) {
                            net_counts[net_id][old_p]--;
                            net_counts[net_id][new_p]++;
                        }
                    }
                }

                for(int i = 0; i < m; ++i) {
                    int span = (net_counts[i][0] > 0) + (net_counts[i][1] > 0) + (net_counts[i][2] > 0);
                    if(span > 1) current_cost++; 
                }
            }

            // FM Local Search Pass
            bool overall_improvement = true;
            
            while(overall_improvement) {
                overall_improvement = false;
                fill(locked.begin(), locked.end(), false);
                history.clear();
                
                fill(head.begin(), head.end(), -1);
                fill(in_bucket.begin(), in_bucket.end(), 0);
                int max_active_gain = 0;

                auto insert_node = [&](int c, int to_p, int g) {
                    int idx = g + OFFSET;
                    if(idx < 0 || idx >= MAX_GAIN) return; 
                    int id = c * 3 + to_p;
                    int curr_head = head[idx];
                    next_node[id] = curr_head;
                    prev_node[id] = -1;
                    if(curr_head != -1) prev_node[curr_head] = id;
                    head[idx] = id;
                    in_bucket[id] = 1;
                    if(idx > max_active_gain) max_active_gain = idx;
                };

                auto remove_node = [&](int c, int to_p, int g) {
                    int id = c * 3 + to_p;
                    if(!in_bucket[id]) return;
                    int idx = g + OFFSET;
                    if(idx < 0 || idx >= MAX_GAIN) return; 
                    int p = prev_node[id];
                    int nxt = next_node[id];
                    if(p != -1) next_node[p] = nxt;
                    else head[idx] = nxt;
                    if(nxt != -1) prev_node[nxt] = p;
                    in_bucket[id] = 0;
                };

                auto update_gain = [&](int cell, int to_part, int delta) {
                    if(locked[cell]) return;
                    int old_g = cell_gain[cell][to_part];
                    int new_g = old_g + delta;
                    remove_node(cell, to_part, old_g);
                    cell_gain[cell][to_part] = new_g;
                    insert_node(cell, to_part, new_g);
                };

                for(int i = 0; i < n; ++i) {
                    int from_p = current_part[i];
                    for(int to_p = 0; to_p < 3; ++to_p) {
                        if(from_p == to_p) continue;
                        int gain = 0;
                        for(int net_id : cells[i].nets) {
                            int S = nets[net_id].cells.size();
                            if (net_counts[net_id][from_p] == S) gain--; 
                            else if (net_counts[net_id][from_p] == 1 && net_counts[net_id][to_p] == S - 1) gain++; 
                        }
                        cell_gain[i][to_p] = gain;
                        insert_node(i, to_p, gain);
                    }
                }

                int step_cost = current_cost;
                int pass_best_cost = current_cost;
                int pass_best_step = -1;

                for(int step = 0; step < n; ++step) {
                    int best_c = -1, best_to_p = -1, best_gain = -1e9;
                    
                    for (int g_idx = max_active_gain; g_idx >= 0; --g_idx) {
                        if (head[g_idx] == -1) {
                            if (g_idx == max_active_gain) max_active_gain--;
                            continue;
                        }
                        bool found = false;
                        int curr_id = head[g_idx];
                        while(curr_id != -1) {
                            int c = curr_id / 3;
                            int to_p = curr_id % 3;
                            int from_p = current_part[c];
                            if (part_sizes[from_p] - 1 < lower_bound || part_sizes[to_p] + 1 > upper_bound) {
                                curr_id = next_node[curr_id];
                                continue;
                            }
                            best_c = c;
                            best_to_p = to_p;
                            best_gain = g_idx - OFFSET;
                            found = true;
                            break;
                        }
                        if (found) break;
                    }

                    if(best_c == -1) break;

                    int from_p = current_part[best_c];
                    for(int t_p = 0; t_p < 3; ++t_p) {
                        if(from_p != t_p) remove_node(best_c, t_p, cell_gain[best_c][t_p]);
                    }
                    locked[best_c] = true;

                    for(int net_id : cells[best_c].nets) {
                        int S = nets[net_id].cells.size();
                        bool critical_before = (net_counts[net_id][0] >= S - 1 || net_counts[net_id][1] >= S - 1 || net_counts[net_id][2] >= S - 1);
                        
                        net_counts[net_id][from_p]--;
                        net_counts[net_id][best_to_p]++;
                        
                        bool critical_after = (net_counts[net_id][0] >= S - 1 || net_counts[net_id][1] >= S - 1 || net_counts[net_id][2] >= S - 1);

                        if (critical_before || critical_after) {
                            for (int neighbor : nets[net_id].cells) {
                                if (neighbor == best_c || locked[neighbor]) continue;
                                int px = current_part[neighbor];
                                
                                for (int k = 0; k < 3; ++k) {
                                    if (k == px) continue;
                                    
                                    int c_px_new = net_counts[net_id][px];
                                    int c_k_new  = net_counts[net_id][k];
                                    int new_g = 0;
                                    if (c_px_new == S) new_g = -1;
                                    else if (c_px_new == 1 && c_k_new == S - 1) new_g = 1;
                                    
                                    int c_px_old = c_px_new + (px == from_p ? 1 : (px == best_to_p ? -1 : 0));
                                    int c_k_old  = c_k_new  + (k == from_p ? 1 : (k == best_to_p ? -1 : 0));
                                    int old_g = 0;
                                    if (c_px_old == S) old_g = -1;
                                    else if (c_px_old == 1 && c_k_old == S - 1) old_g = 1;
                                    
                                    if (new_g != old_g) {
                                        update_gain(neighbor, k, new_g - old_g);
                                    }
                                }
                            }
                        }
                    }
                    
                    part_sizes[from_p]--;
                    part_sizes[best_to_p]++;
                    current_part[best_c] = best_to_p;
                    step_cost -= best_gain;

                    history.push_back({best_c, from_p, best_to_p, best_gain});

                    if(step_cost < pass_best_cost) {
                        pass_best_cost = step_cost;
                        pass_best_step = step;
                    }
                }

                if(pass_best_cost >= current_cost) {
                    for(int i = (int)history.size() - 1; i >= 0; --i) {
                        int c = history[i].cell;
                        int from_p = history[i].from;
                        int to_p = history[i].to;
                        for(int net_id : cells[c].nets) {
                            net_counts[net_id][to_p]--;
                            net_counts[net_id][from_p]++;
                        }
                        part_sizes[to_p]--;
                        part_sizes[from_p]++;
                        current_part[c] = from_p;
                    }
                    break;
                } else {
                    for(int i = (int)history.size() - 1; i > pass_best_step; --i) {
                        int c = history[i].cell;
                        int from_p = history[i].from;
                        int to_p = history[i].to;
                        for(int net_id : cells[c].nets) {
                            net_counts[net_id][to_p]--;
                            net_counts[net_id][from_p]++;
                        }
                        part_sizes[to_p]--;
                        part_sizes[from_p]++;
                        current_part[c] = from_p;
                    }
                    current_cost = pass_best_cost;
                    overall_improvement = true;
                }
            }

            // SA acceptance criteria
            if (current_cost < accepted_cost) {
                accepted_cost = current_cost;
                accepted_part = current_part;
            } else {
                double prob = exp(static_cast<double>(accepted_cost - current_cost) / T);
                if (dist_prob(rng) < prob) {
                    accepted_cost = current_cost;
                    accepted_part = current_part;
                }
            }

            if(current_cost < best_overall_cut) {
                best_overall_cut = current_cost;
                best_overall_part = current_part;
                non_improve_iters = 0; 
            } else {
                non_improve_iters++;
            }

            // Temperature decay
            T *= alpha; 
            if (T < 0.01) T = 100.0; 

            // Restart if trapped
            if (non_improve_iters > 50) {
                force_restart = true;
                non_improve_iters = 0;
                T = 100.0;
            }
        }

        vector<string> g1_cells, g2_cells, g3_cells;
        for(int i = 0; i < n; ++i) {
            if(best_overall_part[i] == 0) g1_cells.push_back(cells[i].name);
            else if(best_overall_part[i] == 1) g2_cells.push_back(cells[i].name);
            else g3_cells.push_back(cells[i].name);
        }

        sort(g1_cells.begin(), g1_cells.end());
        sort(g2_cells.begin(), g2_cells.end());
        sort(g3_cells.begin(), g3_cells.end());

        fout << "Cutsize = " << best_overall_cut << "\n";
        fout << "G1 " << g1_cells.size() << "\n";
        for(size_t i = 0; i < g1_cells.size(); ++i) fout << (i > 0 ? " " : "") << g1_cells[i];
        fout << ";\n";

        fout << "G2 " << g2_cells.size() << "\n";
        for(size_t i = 0; i < g2_cells.size(); ++i) fout << (i > 0 ? " " : "") << g2_cells[i];
        fout << ";\n";

        fout << "G3 " << g3_cells.size() << "\n";
        for(size_t i = 0; i < g3_cells.size(); ++i) fout << (i > 0 ? " " : "") << g3_cells[i];
        fout << ";\n";

        return 0;
    }
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if(argc != 3) {
        cerr << "Usage: ./Lab1 [input file name] [output file name]\n";
        return 1;
    }

    // Fast scan to count nets
    ifstream fin_scan(argv[1]);
    if(!fin_scan.is_open()) {
        cerr << "Error: Cannot open input file.\n";
        return 1;
    }

    int net_count = 0;
    string token;
    while(fin_scan >> token) {
        if(token == "NET") {
            net_count++;
        }
    }
    fin_scan.close();

    cout << "[Pre-Scan] Total Nets detected: " << net_count << endl;

    // Dispatch solver based on net count
    if (net_count < 50000) {
        cout << "[Strategy] Using Pure SA..." << endl;
        return SmallSolver::run(argc, argv);
    } else {
        cout << "[Strategy] Using FM + global SA..." << endl;
        return LargeSolver::run(argc, argv);
    }
}