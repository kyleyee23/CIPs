#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "globaldata.h"
#include "wignerSymbols.h"

using namespace std;

#define pi      3.14159265358979323846  /* pi */
#define ellmax  200
#define TT      1

namespace GlobalData
{
    vector<double> CltildeTT;
    vector<double> ClTT;
    vector<double> ClTTv2;
    vector<double> ClTdT;
    vector<double> Clff;
}

//here, index 1 represents TT
//need to verify this function
//note that these values correspond to l \in (2,ellmax)
vector<double> read_dat(string file_name, int index, bool do_normalizing) {
    ifstream file(file_name.c_str());

    string line;
    vector<double> data;
    int flag = 0;
    while(getline(file, line)){
        if (line.at(0) == '#') {
            getline(file,line);
        }
        vector<double> lineData;
        stringstream lineStream(line);

        double value;
        // Read an integer at a time from the line
        while(lineStream >> value) {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        if (do_normalizing) {
            data.push_back(2*pi*lineData.at(index)/(lineData.at(0) * (lineData.at(0) + 1)));
        }
        else {
            data.push_back(lineData.at(index));
        }
    }
    file.close();
    return data;
}

double F(int l1, int L, int l2) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double return_val = WignerSymbols::wigner3j((double)l1,(double)L,(double)l2,0,0,0);
        return (L*(L + 1) + l2*(l2 + 1) - l1*(l1 + 1))*sqrt(((2*L + 1)*(2*l1 + 1)*(2*l2 + 1))/(16*pi))* return_val;
    }
    else {
        return 0;
    }
}
double fTT(int l1, int L, int l2) {
    return GlobalData::CltildeTT.at(l1 - 2)*F(l2, L, l1) + GlobalData::CltildeTT.at(l2 - 2)*F(l1, L, l2);
}
double gTT(int l1, int L, int l2) {
    return fTT(l1, L, l2)/(2*GlobalData::ClTT.at(l1 - 2)*GlobalData::ClTT.at(l2 - 2));
}
double ATT(int L) {
    double sum = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum += gTT(l1,L,l2)*fTT(l1,L,l2);
        }
    }
    return (L*(L + 1)*(2*L + 1))/sum;
}
double hTT(int l1, int L, int l2) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double return_val = WignerSymbols::wigner3j((double)l1,(double)L,(double)l2,0,0,0);
        return sqrt(((2*L + 1)*(2*l1 + 1)*(2*l2 + 1))/(4/pi))*(GlobalData::ClTdT.at(l1 - 2) + GlobalData::ClTdT.at(l2 - 2))* return_val;
    }
    else {
        return 0;
    }
}

double QTT(int L) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum1 += hTT(l1,L,l2)*gTT(l1,L,l2);
            sum2 += fTT(l1,L,l2)*gTT(l1,L,l2);
        }
    }
    return sum1/sum2;
}

int main() {
    int l = 0;

    GlobalData::CltildeTT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_no_lensing_scalCls.dat", TT, 1);
    GlobalData::ClTT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_lensing_scalCls.dat", TT, 1);
    GlobalData::ClTTv2 = read_dat("../../Dropbox/CIPs MCMC/cmb_files/clXX.dat", TT, 1);
    GlobalData::ClTdT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_X_dX.dat", TT, 1);
    GlobalData::Clff = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_lensing_lenspotentialCls.dat", 5, 0);

    vector<double> result;
    for (int i = 1; i < 15; i++) {
        result.push_back(QTT(i)/pi);
        cout << result.at(i-1) << endl;
    }

    return 0;
}
