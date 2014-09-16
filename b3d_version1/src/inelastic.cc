#ifndef __INELASTIC_CC__
#define __INELASTIC_CC__
#include "b3d.h"

using namespace std;

CB3D *CInelasticList::b3d=NULL;
bool CInelasticList::UseFile = false;
bool CInelasticList::UseInelasticArray = false;

CInelasticList::CInelasticList(){
	string filename;
	ifstream inelasticfile;
	// UseFile = false;
	// UseInelasticArray = true;

	if(b3d!=NULL){
		parmap=&(b3d->parmap);
		GfirstResInfoptr=b3d->reslist->GfirstResInfoptr;
		NResonances = b3d->reslist->NResonances;
		filename = parameter::getS(*parmap,"B3D_INELASTIC_INFO_FILE",string("inelastic.tmp"));

		//create a new NResonances x NResonances array to store lists of elements
		InelasticArray = new list<CInelasticInfo> *[NResonances];
		for(int i = 0; i<NResonances; i++){
			InelasticArray[i] = new list<CInelasticInfo> [NResonances];
		}
		if(UseFile){
			inelasticfile.open(filename.c_str());
			if(inelasticfile){
				cout << "Importing inelastic information from file" << endl;
				ReadInelasticInfo(true);
			}
			else{
				cout << "No inelastic information file, creating..." << endl;
				ReadInelasticInfo(false);
			}
		}
		else{
			cout << "Not Using File" << endl;
			ReadInelasticInfo(false);
		}
	}
}

CInelasticInfo::CInelasticInfo(CResInfo *resinfo_1_in, CResInfo *resinfo_2_in, int type_in){
	//by convention, resinfo_1 will have the lighter of the two particles
	if(resinfo_1_in->mass <= resinfo_2_in->mass){
		resinfo_1 = resinfo_1_in;
		resinfo_2 = resinfo_2_in;
	}
	else{
		resinfo_1 = resinfo_2_in;
		resinfo_2 = resinfo_1_in;
	}
	type = type_in;
	min_mass = resinfo_1->mass + resinfo_2->mass;
	net_q = resinfo_1->charge + resinfo_2->charge;
	net_s = resinfo_1->strange + resinfo_2->strange;
	net_b = resinfo_1->baryon + resinfo_2->baryon;
}

void CInelasticInfo::Print(){
	resinfo_1->Print();
	resinfo_2->Print();
}

void CInelasticList::ReadInelasticInfo(bool FromFile){
	/*
	NOTE: The thermal array indices are named the convention {absolute value of quantities}{sign of quantites}
	and the quantities are named in alphabetical order; baryon number (b), charge (q), strangeness (s).
	The sign function returns 0 if the number is negative, and 1 otherwise.
	*/

	string filler;
	fstream inelasticfile;
	int ires1, ires2, ires3, ires4, netq, netb, nets, pmq=0, pmb=0, pms=0, size, foobar=0, sum = 0;
	double foo = 0;
	CResInfo *resinfoptr_1 = NULL,*resinfoptr_2 = NULL;
	CInelasticInfo *temp, *temp2;
	list<CInelasticInfo>::iterator Th_iter;

	if(FromFile){
		inelasticfile.open(filename.c_str(), fstream::in);
		if(inelasticfile){
			inelasticfile >> filler; //catches #THERMAL tag
			while(sum != 18 && inelasticfile.good() ){
				inelasticfile >> netb >> netq >> nets >> pmb >> pmq >> pms >> size;
				sum = netq + netb + nets + pmq + pmb + pms;
				for(int i=0; i<size; i++){
					inelasticfile >> ires1 >> ires2;
					b3d->reslist->GetResInfoptr(ires1, resinfoptr_1);
					b3d->reslist->GetResInfoptr(ires2, resinfoptr_2);

					CInelasticInfo temp2(resinfoptr_1,resinfoptr_2, 0);
					ThermalArray[netb][netq][nets][pmb][pmq][pms].push_back(temp2);
					foobar++;
				}
			}
			if(UseInelasticArray){
				inelasticfile >> filler;
				if(strcmp(filler.c_str(), "#INELASTIC_ARRAY") == 0){
					while(inelasticfile >> ires1 >> ires2 >> size){
						for(int i=0; i<size; i++){
							inelasticfile >>ires3 >> ires4;
							b3d->reslist->GetResInfoptr(ires3, resinfoptr_1);
							b3d->reslist->GetResInfoptr(ires4, resinfoptr_2);

							CInelasticInfo temp2(resinfoptr_1,resinfoptr_2, 0);
							InelasticArray[ires1][ires2].push_back(temp2);
						}
					}
				}else{
					cout << "Inelastic array needed, but not in inelastic file. Erasing and creating a new one." << endl;
					inelasticfile.close();
					remove(filename.c_str());
					ReadInelasticInfo(false);
					return;
				}
			}
		} else {
			cout << "Error; inelastic file can't be opened" << endl;
			exit(1);
		}
		inelasticfile.close();
	}
	else{
		resinfoptr_1 = GfirstResInfoptr;
			//cout << "Working" << endl;
		while(resinfoptr_1 != NULL){
			resinfoptr_2= resinfoptr_1;
			while(resinfoptr_2 != NULL){
					//create a new CInelasticInfo object from the two current resonance info
				temp = new CInelasticInfo(resinfoptr_1, resinfoptr_2, 0);
				SortedAdd(ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)], *temp);
				foobar++;
				delete temp;
				resinfoptr_2 = resinfoptr_2->nextResInfoptr;
			}
			resinfoptr_1 = resinfoptr_1->nextResInfoptr;
		}
		if(UseInelasticArray){
				//cout << "Creating inelastic array" << endl;
			resinfoptr_1 = GfirstResInfoptr;
			while(resinfoptr_1 != NULL){
				resinfoptr_2 = GfirstResInfoptr;
				while(resinfoptr_2 != NULL){
					cout << resinfoptr_1->ires << "," << resinfoptr_2->ires << endl;
					temp = new CInelasticInfo(resinfoptr_1, resinfoptr_2, 0);
					Th_iter = ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)].begin();
					while(Th_iter != ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)].end()){
						ires1 = Th_iter->resinfo_1->ires;
						ires2 = Th_iter->resinfo_2->ires;
						if(AddToArrayCheck(*(Th_iter->resinfo_1), *(Th_iter->resinfo_2), *resinfoptr_1, *resinfoptr_2)){
							SortedAdd(InelasticArray[min(ires1, ires2)][max(ires1, ires2)], *temp);
							foo++;
								//cout << foo << endl;
							if(foo >= 5000000000){
								cout << "Error, too many resonances are being read in." << endl;
								exit(-1);
							}
						}
						Th_iter++;

					}
					delete temp;
					resinfoptr_2 = resinfoptr_2->nextResInfoptr;
				}
				resinfoptr_1 = resinfoptr_1->nextResInfoptr;
			}
		}
		if(UseFile){
			//cout << "storing to file" << endl;
			inelasticfile.open(filename.c_str(), fstream::out);
			if(inelasticfile){
				//inelasticfile << b3d->run_name << endl;
				inelasticfile << "#THERMAL" << endl;
				for(int i = 0; i<3; i++){
					for(int j = 0; j<5; j++){
						for(int k = 0; k<7; k++){
							for(int pm1 = 0; pm1 <2; pm1++){
								for(int pm2 = 0; pm2<2; pm2++){
									for(int pm3 = 0; pm3<2; pm3++){
										inelasticfile << i <<"\t"<< j << "\t"<< k << "\t"<< pm1 << "\t" << pm2 <<"\t"<< pm3 << "\t" << ThermalArray[i][j][k][pm1][pm2][pm3].size() << endl;
										Th_iter = ThermalArray[i][j][k][pm1][pm2][pm3].begin();
										while(Th_iter != ThermalArray[i][j][k][pm1][pm2][pm3].end()){
											inelasticfile << "\t" << Th_iter->resinfo_1->code << "\t" << Th_iter->resinfo_2->code << endl;
											Th_iter++;
										}
										inelasticfile << endl;
									}
								}
							}
						}
					}
				}
				if(UseInelasticArray){
					inelasticfile << "#INELASTIC_ARRAY" << endl;
					for(int i = 1; i <= b3d->reslist->NResonances; i++){
						for(int j = 1; j <= b3d->reslist->NResonances; j++){
							list<CInelasticInfo> temp = InelasticArray[i][j];
							if(temp.size() > 0){
								inelasticfile << i << "\t" << j << "\t" << temp.size() << endl;
								Th_iter = InelasticArray[i][j].begin();
								while(Th_iter != InelasticArray[i][j].end()){
									inelasticfile << "\t" << Th_iter->resinfo_1->ires << "\t" << Th_iter->resinfo_2->ires << endl;
								}
							}
						}
					}
				}
			}else{
				cout << "Unable to open inelastic file " << filename << endl;
				exit(-1);
			}
			inelasticfile.close();	
		}
	}
	cout << "Stored " << foobar << " inelastic resonances" << endl;
}

// adds CInelasticInfo inelastic_in to list list_in, assuming that list_in is sorted by energy
void CInelasticList::SortedAdd(list<CInelasticInfo> &list_in, CInelasticInfo inelastic_in){
	list<CInelasticInfo>::iterator CLiterator = list_in.begin();

	while(CLiterator->min_mass < inelastic_in.min_mass && CLiterator != list_in.end() ){
		CLiterator++;
	}
	list_in.insert(CLiterator, inelastic_in);
}

/*
Method to check whether or not a given inelastic scattering is valid. Allows for more complex checks.
Assumes that res1 and res2 are the incoming particles, while res3 and res4 are outgoing particles.
*/
bool CInelasticList::AddToArrayCheck(CResInfo res1, CResInfo res2, CResInfo res3, CResInfo res4){
	bool add = true;
	//check for charge, baryon number, and strangness conservation
	if(res1.charge + res2.charge - res3.charge - res4.charge != 0){
		add = false;
	}
	if(res1.baryon + res2.baryon - res3.baryon - res4.baryon != 0){
		add = false;
	}
	if(res1.strange + res2.strange - res3.strange - res4.strange != 0){
		add = false;
	}

	//check that the incoming and outgoing resonances aren't the same
	if((res1.ires == res3.ires && res2.ires == res4.ires)||(res1.ires == res4.ires && res2.ires == res3.ires)){
		add = false;
	}

	return add;
}

#endif