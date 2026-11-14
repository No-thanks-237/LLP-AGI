/* drv.cpp -- batch driver for the AUTHOR'S Light_scalar_decay CP-even-H module.
 *
 * PATCH-2 (decay module): the original T1Aa.cpp main() (not shipped in
 * inputs/decay_source) is replaced by this batch driver.  NO physics is
 * changed: for every (mH, tanb) the model is initialised with exactly the
 * sequence prescribed by manual.pdf / prompt_2HDM.md section 3
 *   set_model(1); set_mH(mH); set_mA(600); set_mHC(600);
 *   set_lambdav(0); set_tanb(tanb); set_cosba(1/tanb);
 * and the public API of Model.h is read.
 *
 * stdin : one point per line "<mH> <tanb>"
 * stdout: one line per point, 13 whitespace separated numbers (%.10e):
 *   mH tanb GammaTot Gamma_gaga Gamma_ee Gamma_mumu Gamma_tautau
 *           Gamma_pipi Gamma_KK Gamma_4pi Gamma_cc Gamma_ss Gamma_gg
 *   (all widths in GeV)
 */
#include <cstdio>
#include <iostream>
#include <string>

#include "Model.h"

int main()
{
    BSModel model;
    double mH_in, tb_in;
    char line[512];

    while (std::cin.getline(line, sizeof(line)))
    {
        std::string s(line);
        if (s.empty() || s[0] == '#') continue;
        if (std::sscanf(s.c_str(), "%lf %lf", &mH_in, &tb_in) != 2) continue;
        if (!(mH_in > 0.0) || !(tb_in > 0.0)) continue;

        model.set_model(1);
        model.set_mH(mH_in);
        model.set_mA(600.0);
        model.set_mHC(600.0);
        model.set_lambdav(0.0);
        model.set_tanb(tb_in);
        model.set_cosba(1.0 / tb_in);

        double Gtot   = model.DecayWidth_HTotal();
        double Ggaga  = model.DecayWidth_Hgaga();
        double Gee    = model.DecayWidth_Hee();
        double Gmumu  = model.DecayWidth_Hmumu();
        double Gtata  = model.DecayWidth_Htautau();
        double Gpipi  = model.DecayWidth_HPiPi();
        double GKK    = model.DecayWidth_HKK();
        double G4pi   = model.DecayWidth_HPiPiPiPi();
        double Gcc    = model.DecayWidth_Hcc();
        double Gss    = model.DecayWidth_Hss();
        double Ggg    = model.DecayWidth_HGluon();

        std::printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                    mH_in, tb_in, Gtot, Ggaga, Gee, Gmumu, Gtata,
                    Gpipi, GKK, G4pi, Gcc, Gss, Ggg);
        std::fflush(stdout);
    }
    return 0;
}
