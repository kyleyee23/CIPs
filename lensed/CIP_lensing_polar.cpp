#include "CIP_lensing_polar.h"

int main() {
    GlobalData::CltildeTT = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", TT, 1);
    GlobalData::ClTT = read_dat("../cmb_files/Cl_lensing_scalCls.dat", TT, 1);
    GlobalData::ClTdT = read_dat("../cmb_files/Cl_X_dX.dat", TT, 1);
    GlobalData::Clff = read_dat("../cmb_files/Cl_lensing_lenspotentialCls.dat", 5, 0);

    GlobalData::CltildeEE = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", EE, 1);
    GlobalData::ClEE = read_dat("../cmb_files/Cl_lensing_scalCls.dat", EE, 1);
    GlobalData::ClEdE = read_dat("../cmb_files/Cl_X_dX.dat", 4, 1);

    GlobalData::CltildeTE = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", 3, 1);
    GlobalData::ClTE = read_dat("../cmb_files/Cl_lensing_scalCls.dat", 4, 1);
    GlobalData::ClTdE = read_dat("../cmb_files/Cl_X_dX.dat", 2, 1);
    GlobalData::ClEdT = read_dat("../cmb_files/Cl_X_dX.dat", 3, 1);

    GlobalData::CltildeTB = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", 4, 1);
    GlobalData::ClTB = read_dat("../cmb_files/Cl_lensing_scalCls.dat", 3, 1);

    GlobalData::CltildeBB = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", 3, 1);
    GlobalData::ClBB = read_dat("../cmb_files/Cl_lensing_scalCls.dat", 4, 1);

    GlobalData::CltildeBE = read_dat("../cmb_files/Cl_no_lensing_scalCls.dat", 3, 1);
    GlobalData::ClBE = read_dat("../cmb_files/Cl_lensing_scalCls.dat", 4, 1);

    double thetaFWHM=7.1/60;
    double sigmab= 0.00741173*thetaFWHM;
    double fsurv = 0.65;
    double NET = 62;
    double tobs = 1.2*3E7;

    for(int l1=2;l1<ellmax; l1++){
      GlobalData::ClTT.at(l1 - 2) =  (GlobalData::ClTT.at(l1 - 2))*exp(-l1*(l1+1)*pow(sigmab,2))+4*pi*fsurv*pow(NET,2)/tobs;
      GlobalData::ClEE.at(l1 - 2) =  (GlobalData::ClEE.at(l1 - 2))*exp(-l1*(l1+1)*pow(sigmab,2))+4*pi*fsurv*pow(NET,2)/tobs;
      GlobalData::ClTE.at(l1 - 2) =  (GlobalData::ClTE.at(l1 - 2))*exp(-l1*(l1+1)*pow(sigmab,2))+4*pi*fsurv*pow(NET,2)/tobs;
    }

    vector<double> result;
    for (int i = 1; i <= 20; i++) {
       result.push_back(QEE(i));
       cout << i ;
       cout << ',';
       cout << result.at(i-1) << endl;
    }

    //int L=500;


    //cout << QTT(L)/pi << endl;

    //int ellmax=2000;

    //double sum = 0.0;

    //for (int l1 = 2; l1 <= ellmax; l1++) {
    // for (int l2 = 2; l2 <= ellmax; l2++) {
    //sum +=  sqrt(((2*(double)L + 1)*(2*(double)l1 + 1)*(2*(double)l2 + 1))/(16*pi));
    //}
    //}

    //cout << sum << endl;

    return 0;
}
