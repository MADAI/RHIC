#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <boost/functional/hash.hpp> // needed for hashing pairs

#include "b3d.h"
#include "hist.h"
#include "pow.h"

int sign(double x) {
    return (x > 0) - (x < 0);
}

/// Calculate the balance functions.
void CB3D::CalcBalance() {
    std::cout << "*** Calculating balance functions for run " << run_name << "..." << std::endl;
    NSAMPLE = parameter::getI(parmap, "B3D_NSAMPLE", 1);
    const int neventsmax = parameter::getI(parmap, "B3D_NEVENTSMAX", 10);

    const string dirname = "analysis/" + run_name + "/" + qualifier;
    const string command = "mkdir -p " + dirname;
    system(command.c_str());
    const string anal_output_filename = dirname + "/results_balance.dat";
    auto anal_output = fopen(anal_output_filename.c_str(), "w");

    parameter::set(parmap, string("B3D_RESONANCES_DECAYS_FILE"), string("progdata/madai/resinfo/decays_pdg_weak.dat"));
    parameter::set(parmap, string("B3D_RESONANCES_INFO_FILE"), string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
    reslist = new CResList();
    reslist->ReadResInfo();

    const double DELRAPIDITY = 2.0 * parameter::getD(parmap, "B3D_ETAMAX", 1.0);
    COLLISIONS = false;
    oscarfile = NULL;
    //^ We need to set this, so the header gets read. Otherwise `feof(oscarfile) == true`.

    /// Balance functions to calculate.
    // (pi+, pi-), (K+, K-), (p, pbar)
    vector<pair<int, int>> balance_pairs = { {211, -211}, {321, -321}, {2212, -2212} };
    /// We need to calculate the histograms for (a,-b), (-a,b) and (-a,-b) too.
    const size_t old_size = balance_pairs.size();
    for (size_t i = 0; i < old_size; i++) {
        const auto a = balance_pairs[i].first;
        const auto b = balance_pairs[i].second;
        balance_pairs.emplace_back(make_pair(-a, b));
        balance_pairs.emplace_back(make_pair(a, -b));
        balance_pairs.emplace_back(make_pair(-a, -b));
    }
    /// Relevant species (ignore anti particles).
    unordered_set<int> balance_species;
    for (const auto& pair : balance_pairs) {
        balance_species.insert(abs(pair.first));
        balance_species.insert(abs(pair.second));
    }
    /// A rapidity histogram for each balance function.
    unordered_map<pair<int, int>, Histogram<double>, boost::hash<pair<int, int>>> histograms;
    // TODO: figure out what the following constants should be.
    const size_t nbins = 80;
    const double min_y = -1;
    const double max_y = 1;
    const double max_dy = max_y - min_y;
    const double min_dy = 0;
    const double min_deta = 0;
    const double max_deta = 7;
    for (const auto& pair : balance_pairs) {
        histograms.emplace(make_pair(pair, Histogram<double>(nbins, min_dy, max_dy)));
    }
    /// A pseudo-rapidity histogram for the charge-only balance function.
    auto charge_histogram = Histogram<double>(nbins, min_deta, max_deta);

    /// Number for each species.
    /// (Needed for denominator of balance function.)
    unordered_map<int, int64_t> total_number(balance_species.size());
    for (const auto& PID : balance_species) {
        total_number[PID] = 0;
    }
    /// Total number including all species.
    int64_t total_number_all = 0;

    int nevents = 0;
    do {
        KillAllParts();
        const int nparts = ReadOSCAR(nevents + 1);
        std::cout << "Processing " << nparts << " particles..." << std::endl;
        if (nparts > 0) {
            nevents += 1;
            //printf("before, nparts=%d =? %d\n",nparts,int(FinalPartMap.size()));
            PerformAllActions(); // Decays unstable particles
            //printf("after decays, nparts=%d\n",int(FinalPartMap.size()));

            // We are only interested in particles occuring in the balance function
            // and with y_min <= y <= y_max.
            vector<CPart> relevant_particles; // TODO: only remember weight, charge, rapidity and momentum
            for (const auto& ppos : FinalPartMap) {
                const auto part = ppos.second;
                const int PID = part->resinfo->code;
                const int charge = part->resinfo->charge;
                const int weight = part->weight;
                assert(weight == 1);
                const double y = part->y;
                // Skip particles we do not care about.
                if (balance_species.count(abs(PID)) == 0 || y < min_y || y > max_y)
                    continue;
                // TODO: p_T cuts and efficiency
                relevant_particles.push_back(*part);
                total_number[abs(PID)] += weight * abs(charge);
                total_number_all += weight * abs(charge);
            }
            std::cout << "Relevant particles: " << relevant_particles.size() << std::endl;

            // Iterate over all pairs to calculate dy and fill histograms.
            for (size_t i = 0; i < relevant_particles.size(); i++) for (size_t j = 0; j < i; j++) {
                const auto& p_i = relevant_particles[i];
                const auto& p_j = relevant_particles[j];
                const double dy = fabs(p_i.y - p_j.y);
                const auto mypair = make_pair(p_i.resinfo->code, p_j.resinfo->code);
                const auto mypair_reversed = make_pair(mypair.second, mypair.first);

                // For B(ab): Calculate N_a N_b.
                auto search = histograms.find(mypair);
                if (search != histograms.end()) {
                    const int charge = p_j.resinfo->charge;
                    const int weight = p_j.weight;
                    search->second.add(dy, weight * abs(charge));
                }
                search = histograms.find(mypair_reversed);
                if (search != histograms.end()) {
                    const int charge = p_i.resinfo->charge;
                    const int weight = p_i.weight;
                    search->second.add(dy, weight * abs(charge));
                }

                // For B(+-).
                if (sign(p_i.resinfo->charge) == -sign(p_j.resinfo->charge)) {
                    double p_i_mag = sqrt(square(p_i.p[1]) + square(p_i.p[2]) + square(p_i.p[3]));
                    double p_j_mag = sqrt(square(p_j.p[1]) + square(p_j.p[2]) + square(p_j.p[3]));
                    // Avoid dividing by zero (which would result in nans).
                    if (p_i_mag == 0)
                        p_i_mag = 1;
                    if (p_j_mag == 0)
                        p_j_mag = 1;
                    const double dpseudorapidity = fabs( atanh(p_i.p[3] / p_i_mag) - atanh(p_j.p[3] / p_j_mag ) );
                    //assert(dpseudorapidity == dpseudorapidity);
                    // B_+- ~ N_+ N_- - N_+ N_+ - N_- N_- + N_- N_+
                    charge_histogram.add(dpseudorapidity,
                        -p_i.resinfo->charge * p_i.weight * p_j.resinfo->charge * p_j.weight);
                }
            }
        }
    } while (!feof(oscarfile) && nevents < neventsmax);
    if (nevents != neventsmax) {
        printf("EVENT SHORTAGE? Got %i, expected %i.\n", nevents, neventsmax);
        std::cout << oscarfile << " vs. " << feof(oscarfile) << std::endl;
    }
    fclose(oscarfile);
    oscarfile = NULL;

    // Calculate balance functions.
    // B_ab(y) = (N_a(0) - N_-a(0))(N_b(y) - N_-b(y)) / (N_b(y) + N_-b(y))
    for (const auto& mypair : balance_pairs) {
        const auto& histogram = histograms.find(mypair)->second;
        vector<int64_t> B(nbins, 0);
        const int a = mypair.first;
        const int b = mypair.second;
        // We need to calculate N_a N_b - N_a N_-b - N_-a N_b + N_-a N_-b.
        const vector<pair<int, int>> pairs = { {a, b}, {a, -b}, {-a, b}, {-a, -b} };
        for (size_t i = 0; i < B.size(); i++) {
            B[i] += histogram.histogram[i];
            B[i] -= histograms.find(pairs[1])->second.histogram[i];
            B[i] -= histograms.find(pairs[2])->second.histogram[i];
            B[i] += histograms.find(pairs[3])->second.histogram[i];
        }
        fprintf(anal_output, "B(%i, %i)\n", a, b);
        fprintf(anal_output, "y  B\n");
        for (size_t i = 0; i < B.size(); i++) {
            const double y = histogram.get_value(i);
            //^ Does not matter which histogram we use.
            fprintf(anal_output, "%f  %lli\n", y, B[i]);
        }
        fprintf(anal_output, "total number of %i\n%lli\n", b, total_number[abs(b)]);
        fprintf(anal_output, "\n");
    }

    vector<int64_t> B_charge(nbins);
    fprintf(anal_output, "B(+, -)\n");
    for (size_t i = 0; i < B_charge.size(); i++) {
        B_charge[i] = charge_histogram.histogram[i];
        const double pseudorapidity = charge_histogram.get_value(i);
        fprintf(anal_output, "%f  %lli\n", pseudorapidity, B_charge[i]);
    }
    fprintf(anal_output, "total number\n%lli\n", total_number_all);
    std::cout << "B(+, -): " << (double) charge_histogram.sum()/charge_histogram.nbins << " counts per bin." << std::endl;

    delete reslist;
    fclose(anal_output);
}
