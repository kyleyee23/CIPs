#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "wignerSymbols.h"
using namespace std;

#define pi      3.14159265358979323846  /* pi */
#define ellmax  200

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
double fTT(int l1, int L, int l2, vector<double> CltildeTT) {
    return CltildeTT.at(l1 - 2)*F(l2, L, l1) + CltildeTT.at(l2 - 2)*F(l1, L, l2);
}
double gTT(int l1, int L, int l2, vector<double> ClTT, vector<double> CltildeTT) {
    return fTT(l1, L, l2,CltildeTT)/(2*ClTT.at(l1 - 2)*ClTT.at(l2 - 2));
}
double ATT(int L, vector<double> ClTT, vector<double> CltildeTT) {
    double sum = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum += gTT(l1,L,l2, ClTT, CltildeTT)*fTT(l1,L,l2, CltildeTT);
        }
    }
    return (L*(L + 1)*(2*L + 1))/sum;
}
double hTT(int l1, int L, int l2, vector<double> ClTdT) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double return_val = WignerSymbols::wigner3j((double)l1,(double)L,(double)l2,0,0,0);
        return sqrt(((2*L + 1)*(2*l1 + 1)*(2*l2 + 1))/(4/pi))*(ClTdT.at(l1 - 2) + ClTdT.at(l2 - 2))* return_val;
    }
    else {
        return 0;
    }
}

double QTT(int L, vector<double> ClTdT, vector<double> ClTT, vector<double> CltildeTT) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum1 += hTT(l1,L,l2, ClTdT)*gTT(l1,L,l2,ClTT,CltildeTT);
            sum2 += fTT(l1,L,l2, CltildeTT)*gTT(l1,L,l2,ClTT,CltildeTT);
        }
    }
    return sum1/sum2;
}

int main() {
    int l = 0;
    int TT = 1;
    vector<double> CltildeTT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_no_lensing_scalCls.dat", TT, 1);
    vector<double> ClTT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_lensing_scalCls.dat", TT, 1);
    vector<double> ClTTv2 = read_dat("../../Dropbox/CIPs MCMC/cmb_files/clXX.dat", TT, 1);
    vector<double> ClTdT = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_X_dX.dat", TT, 1);
    vector<double> Clff = read_dat("../../Dropbox/CIPs MCMC/cmb_files/Cl_lensing_lenspotentialCls.dat", 5, 0);

    vector<double> result;
    for (int i = 1; i < 15; i++) {
        result.push_back(QTT(i, ClTdT, ClTT, CltildeTT)/pi);
        cout << result.at(i-1) << endl;
    }

    return 0;
}