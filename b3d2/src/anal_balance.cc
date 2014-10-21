#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "b3d.h"
#include "hist.h"

/// Represent all balance functions.
/// { (PID, PID): B(0, eta), ... }
//typedef unordered_map<pair<int, int>, Histogram<double>> BalanceMap;

/// Calculate the balance functions.
void CB3D::CalcBalance() {
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

    /// Balance functions to calculate.
    vector<pair<int, int>> balance_pairs = { {211, -211}, {321, -321} };
    /// Relevant species.
    unordered_set<int> balance_species;
    for (const auto& pair : balance_pairs) {
        balance_species.insert(pair.first);
        balance_species.insert(-pair.first);
        balance_species.insert(pair.second);
        balance_species.insert(-pair.second);
    }
    /// A histogram for each species
    unordered_map<int, Histogram<double>> histograms;
    histograms.reserve(balance_species.size());
    // TODO: figure out what the following constants should be.
    const size_t nbins = 100;
    const double min_eta = -2;
    const double max_eta = 2;
    for (const auto& PID : balance_species) {
        histograms.emplace(make_pair(PID, Histogram<double>(nbins, min_eta, max_eta)));
    }

	int nevents = 0;
	do {
	    KillAllParts();
		const int nparts = ReadOSCAR(nevents+1);
		if (nparts > 0){
			nevents += 1;
			//printf("before, nparts=%d =? %d\n",nparts,int(FinalPartMap.size()));
			PerformAllActions(); // Decays unstable particles
			//printf("after decays, nparts=%d\n",int(FinalPartMap.size()));

            // We need a histogram in eta for each species and antispecies.
			for (const auto& ppos : FinalPartMap) {
				const auto part = ppos.second;
				const int PID = part->resinfo->code;
                const auto search = histograms.find(PID);
                if (search != histograms.end()) {
                    const double eta = part->eta;
                    const int weight = part->weight;
                    const int charge = abs(part->resinfo->charge);
                    search->second.add(eta, weight*charge);
                }
			}
		}
	} while (!feof(oscarfile) && nevents < neventsmax);
	if (nevents != neventsmax) {
		printf("EVENT SHORTAGE?\n");
	}
	fclose(oscarfile);
	oscarfile = NULL;

    // Calculate balance function.
    // B_ab(eta) = (N_a(0) - N_-a(0))(N_b(eta) - N_-b(eta)) / (N_b(eta) + N_-b(eta))
    for (const auto& mypair : balance_pairs) {
        vector<double> B(nbins);
        const int a = mypair.first;
        const int b = mypair.second;
        for (size_t i = 0; i < B.size(); i++) {
            const int Na_pos = histograms.at(a).get_count(0);
            const int Na_neg = histograms.at(-a).get_count(0);
            const int Nb_pos = histograms.at(b).histogram[i];
            const int Nb_neg = histograms.at(-b).histogram[i];
            B[i] = (Na_pos - Na_neg) * (Nb_pos - Nb_neg) / (Nb_pos + Nb_neg);
        }
        fprintf(anal_output, "B(%i, %i)\n", a, b);
        fprintf(anal_output, "eta  B\n");
        for (size_t i = 0; i < B.size(); i++) {
            const double eta = min_eta + (max_eta - min_eta) / nbins * (i + 0.5);
            fprintf(anal_output, "%f  %f", eta, B[i]);
        }
    }

	delete reslist;
	fclose(anal_output);
}
