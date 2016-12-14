#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include <cmath>

#include "globaldata.h"
#include "wignerSymbols.h"

using namespace std;

#define pi      3.14159265358979323846  /* pi */
#define ellmax  2500
#define TT      1
#define gamma   1
#define EE      2

namespace GlobalData
{
  vector<double> CltildeTT;
  vector<double> ClTT;
  vector<double> ClTdT;
  vector<double> Clff;

  vector<double> CltildeEE;
  vector<double> ClEE;
  vector<double> ClEdE;

  vector<double> CltildeTE;
  vector<double> ClTE;

  //index 1
  vector<double> ClTdE;

  //index 2
  vector<double> ClEdT;

}


vector<double> read_dat(string file_name, int index, bool do_normalizing) {
    ifstream file(file_name.c_str());

    string line;
    vector<double> data;
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
    if ((abs(l1 - L) <= l2) && (l2 <= l1 + L)) {
      double dL = (double)L;
      double dl1 = (double)l1;
      double dl2 = (double)l2;
      double return_val = WignerSymbols::wigner3j(dl1,dL,dl2,0,0,0);
        return (dL*(dL + 1) + dl2*(dl2 + 1) - dl1*(dl1 + 1))*sqrt(((2*dL + 1)*(2*dl1 + 1)*(2*dl2 + 1))/(16*pi))* return_val;
    }
    else {
        return 0;
    }
}

double fTT(int l1, int L, int l2) {
    return GlobalData::CltildeTT.at(l1 - 2)*F(l2, L, l1) + GlobalData::CltildeTT.at(l2 - 2)*F(l1, L, l2);
}

double fEE(int l1, int L, int l2) {
    //O&H had an even written next to this equation; does this mean it only works for even values of l1,l2?
    //also, subscripts on CL?

    if ((l1 + L + l2)%2 == 1) {
        return 0;
    }

    return GlobalData::CltildeEE.at(l1-2)*F(l2, L, l1) + GlobalData::CltildeEE.at(l2-2)*F(l1,L,l2);
}

double fTE(int l1, int L, int l2) {
    //O&H had an even written next to this equation; does this mean it only works for even values of l1,l2?
    //also, subscripts on CL?

    if ((l1 + L + l2)%2 == 1) {
        return 0;
    }

    return GlobalData::CltildeTE.at(l1-2)*F(l2, L, l1) + GlobalData::CltildeTE.at(l2-2)*F(l1,L,l2);
}

double gTT(int l1, int L, int l2) {
    return fTT(l1, L, l2)/(2*GlobalData::ClTT.at(l1 - 2)*GlobalData::ClTT.at(l2 - 2));
}

double gEE(int l1, int L, int l2) {
    return fEE(l1, L, l2)/(2*GlobalData::ClEE.at(l1 - 2)*GlobalData::ClEE.at(l2 - 2));
}

double gTE(int l1, int L, int l2) {
    //replace with complex conjugate

    double x = GlobalData::ClTE.at(l1-2);
    double y = GlobalData::ClTE.at(l2-2);
    double result = x*y;
    double num = (GlobalData::ClTT.at(l2-2) * GlobalData::ClEE.at(l2-2)*fTE(l1, L, l2)) - (pow(-1,(L + l1 + l2))*GlobalData::ClTE.at(l1-2)*GlobalData::ClTE.at(l2-2)*fTE(l2,L,l1));
    return num/(GlobalData::ClTT.at(l1-2) * GlobalData::ClTT.at(l2-2) * GlobalData::ClEE.at(l1-2) * GlobalData::ClEE.at(l2-2) - pow(result,2));
}

double ATT(int L) {
    double sum = 0.0;
    double dL = (double)L;

    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum += gTT(l1,L,l2)*fTT(l1,L,l2);
        }
    }
    return (dL*(dL + 1)*(2*dL + 1))/sum;
}

double NTT(int L) {

  double dL = (double)L;

    double sum = 0.0;

    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
	  sum += gTT(l1,L,l2)*((GlobalData::ClTT.at(l1 - 2))*(GlobalData::ClTT.at(l2 - 2))*gTT(l1,L,l2)\
			       +pow(-1,L+l1+l2)*(GlobalData::ClTT.at(l2 - 2))*(GlobalData::ClTT.at(l1 - 2))*gTT(l2,L,l1));
        }
    }

    return ATT(L)*ATT(L)/(dL*(dL+1)*(2*dL+1))*sum;

}

double hTT(int l1, int L, int l2) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double dL=(double)L;
        double dl1=(double)l1;
        double dl2=(double)l2;
        double return_val = WignerSymbols::wigner3j(dl1,dL,dl2,0,0,0);
        return sqrt(((2*dL + 1)*(2*dl1 + 1)*(2*dl2 + 1))/(4*pi))*(GlobalData::ClTdT.at(l1 - 2) + GlobalData::ClTdT.at(l2 - 2))* return_val;
    }
    else {
        return 0;
    }
}

double hEE(int l1, int L, int l2) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double dL=(double)L;
        double dl1=(double)l1;
        double dl2=(double)l2;
        double return_val = WignerSymbols::wigner3j(dl1,dL,dl2,0,0,0);
        return sqrt(((2*dL + 1)*(2*dl1 + 1)*(2*dl2 + 1))/(4*pi))*(GlobalData::ClEdE.at(l1 - 2) + GlobalData::ClEdE.at(l2 - 2))* return_val;
    }
    else {
        return 0;
    }
}

double hTE(int l1, int L, int l2) {
    if ((abs(l1 - l2) <= L) && (L <= l1 + l2)) {
        double dL=(double)L;
        double dl1=(double)l1;
        double dl2=(double)l2;
        double return_val = WignerSymbols::wigner3j(dl1,dL,dl2,0,0,0);
        return sqrt(((2*dL + 1)*(2*dl1 + 1)*(2*dl2 + 1))/(4*pi))*(GlobalData::ClTdE.at(l1 - 2) + GlobalData::ClEdT.at(l2 - 2))* return_val;
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

double QEE(int L) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum1 += hEE(l1,L,l2)*gEE(l1,L,l2);
            sum2 += fEE(l1,L,l2)*gEE(l1,L,l2);
        }
    }
    return sum1/sum2;
}

double QTE(int L) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int l1 = 2; l1 <= ellmax; l1++) {
        for (int l2 = 2; l2 <= ellmax; l2++) {
            sum1 += hTE(l1,L,l2)*gTE(l1,L,l2);
            sum2 += fTE(l1,L,l2)*gTE(l1,L,l2);
        }
    }
    return sum1/sum2;
}

double kappaTT(int l1, int l2) {
    int sum = 0;

    for (int L = 1; L <=100; L++) {
        sum += (2*L + 1)/(L*L) * hTT(l1,L,l2)*hTT(l1,L,l2);
    }
    return gamma*sum;
}

double KTT(int L) {
    int sum = 0;
    double dL=(double)L;
    for (int l1 = 2; l1 < ellmax; l1++) {
        for (int l2 = 2; l2 < ellmax; l2++) {
            sum += gTT(l1,L,l2)*gTT(l1,L,l2)*kappaTT(l1, l2);
        }
    }
    return ATT(L)*ATT(L)/(dL*(dL+1)*(2*dL+1)) * sum;
}
