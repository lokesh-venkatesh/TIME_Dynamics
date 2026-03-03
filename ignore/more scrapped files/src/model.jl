
module CancerModel

using DifferentialEquations

export tumor_ode!

function tumor_ode!(du, u, p, t)
    # State Variables
    S, Sr, C, Cr, M1, M2, Th1, Th2, Tc, Treg, IL10, IFNg, IL2 = u
    
    # Parameters (unpacking from Dict for readability)
    # In production, use a NamedTuple or ComponentArray for speed
    gs = p["gamma_s"]; ms = p["m_s"]; p1 = p["p1"]; p2 = p["p2"]; ds = p["delta_s"]
    mus = p["mu_s"]; k1 = p["k1"]; tck = p["tck"]; ktc1 = p["ktc1"]
    musr = p["mu_sr"]; k2 = p["k2"]; ktc2 = p["ktc2"]
    gc = p["gamma_c"]; mc = p["m_c"]; K = p["Ktumor"]; r1 = p["r1"]; dc = p["delta_c"]
    muc1 = p["mu_c1"]; k3 = p["k3"]; muc2 = p["mu_c2"]; k4 = p["k4"]; ktc3 = p["ktc3"]
    r2 = p["r2"]; dcr = p["delta_cr"]; muc1r = p["mu_c1_r"]; k5 = p["k5"]; muc2r = p["mu_c2_r"]
    k6 = p["k6"]; ktc4 = p["ktc4"]
    gm1 = p["gamma_m1"]; lm1 = p["lambda_m1"]; dm1 = p["delta_m1"]; mm1ck2 = p["mu_m1ck2"]; k7 = p["k7"]
    gm2 = p["gamma_m2"]; lm2 = p["lambda_m2"]; dm2 = p["delta_m2"]; mm2ck1 = p["mu_m2ck1"]; k10 = p["k10"]
    gth1 = p["gamma_th1"]; lth1 = p["lambda_th1"]; dth1 = p["delta_th1"]; mth1ck1 = p["mu_th1ck1"]; k8 = p["k8"]
    mth1ck3 = p["mu_th1ck3"]; k9 = p["k9"]
    gth2 = p["gamma_th2"]; lth2 = p["lambda_th2"]; dth2 = p["delta_th2"]
    gtc = p["gamma_tc"]; ltc1 = p["lambda_tc1"]; ltc4 = p["lambda_tc4"]; mtcs = p["mu_tcs"]
    ltc2 = p["lambda_tc2"]; dtc = p["delta_tc"]; mtctreg = p["mu_tctreg"]; ltc3 = p["lambda_tc3"]
    gtreg = p["gamma_treg"]; ltreg2 = p["lambda_treg2"]; dtreg = p["delta_treg"]
    mtregck1 = p["mu_tregck1"]; k11 = p["k11"]
    bm2 = p["beta_m2"]; dck1 = p["delta_ck1"]; btreg = p["beta_treg"]; bth2 = p["beta_th2"]
    bth1ck2 = p["beta_th1ck2"]; btc = p["beta_tc"]; dck2 = p["delta_ck2"]
    bth1ck3 = p["beta_th1ck3"]; dck3 = p["delta_ck3"]

    # 1. Stem Cells (S)
    term_renewal = (gs * (1 - ms) * (1 - p1 - p2)) * S
    term_diff_death = (ds + (p2 * gs) + (gs * ms * p1 / 2.0)) * S
    term_immune_s = (mus * S * IFNg) / (k1 + IFNg)
    term_tc_s = (tck * S * Tc) / (ktc1 + Tc)
    du[1] = term_renewal - term_diff_death - term_immune_s - term_tc_s

    # 2. Resistant Stem Cells (Sr)
    term_sr_growth = (gs * (1 - p1 - p2) - (ds + p2 * gs)) * Sr
    term_sr_mut = ms * gs * (1 - p1/2.0 - p2) * S
    term_immune_sr = (musr * Sr * IFNg) / (k2 + IFNg)
    term_tc_sr = (tck * Sr * Tc) / (ktc2 + Tc)
    du[2] = term_sr_growth + term_sr_mut - term_immune_sr - term_tc_sr

    # 3. Cancer Cells (C)
    gompertz_c = gc * (1 - mc) * log(0.5 * K / (C + r1)) * C
    term_c_diff = gs * (p1 + p2) * S
    term_c_death = dc * C + mc * gc * C
    term_immune_c_il10 = (muc1 * C * IL10) / (IL10 + k3)
    term_immune_c_ifng = (muc2 * C * IFNg) / (IFNg + k4)
    term_tc_c = (tck * C * Tc) / (ktc3 + Tc)
    du[3] = gompertz_c + term_c_diff - term_c_death + term_immune_c_il10 - term_immune_c_ifng - term_tc_c

    # 4. Resistant Cancer Cells (Cr)
    gompertz_cr = gc * Cr * log(0.5 * K / (Cr + r2))
    term_cr_diff = gs * Sr * (p1 + p2) + mc * gc * C
    term_cr_death = dcr * Cr
    term_immune_cr_il10 = (muc1r * Cr * IL10) / (IL10 + k5)
    term_immune_cr_ifng = (muc2r * Cr * IFNg) / (IFNgamma + k6)
    term_tc_cr = (tck * Cr * Tc) / (ktc4 + Tc)
    du[4] = gompertz_cr + term_cr_diff - term_cr_death + term_immune_cr_il10 - term_immune_cr_ifng - term_tc_cr

    # 5. M1 Macrophages
    du[5] = gm1 * M1 * (C + Cr) / (M1 + lm1) - dm1 * M1 + (mm1ck2 * M1 * IFNg) / (IFNg + k7)

    # 6. M2 Macrophages
    du[6] = gm2 * M2 * (C + Cr) / (M2 + lm2) - dm2 * M2 + (mm2Ck1 * M2 * IL10) / (IL10 + k10)

    # 7. TH1 Cells
    du[7] = gth1 * (TH1 * M1) / (lth1 + TH1) - dth1 * TH1 - (mth1ck1 * IL10 * TH1) / (IL10 + k8) + (mth1ck3 * IL2 * TH1) / (IL2 + k9)

    # 8. TH2 Cells
    du[8] = gth2 * (TH2 * M2) / (lth2 + TH2) - dth2 * TH2

    # 9. Cytotoxic T Cells (Tc)
    term_tc_act1 = gtc * Tc * (C + Cr) / (Tc + ltc1)
    term_tc_act2 = gtc * Tc * TH1 / (Tc + ltc4)
    term_tc_inh1 = (mtcs * Tc * (S + Sr)) / (Tc + ltc2)
    term_tc_inh2 = (mTcTr_eg * Tc * Treg) / (ltc3 + Treg)
    du[9] = term_tc_act1 + term_tc_act2 - term_tc_inh1 - dtc * Tc - term_tc_inh2

    # 10. Regulatory T Cells (Treg)
    du[10] = gtreg * (Treg * M2) / (Treg + ltreg2) - dtreg * Treg + (mtregck1 * IL10 * Treg) / (Treg + k11)

    # 11. IL-10 Cytokine
    du[11] = bm2 * M2 - dck1 * IL10 + btreg * Treg + bth2 * TH2

    # 12. IFN-gamma
    du[12] = bth1ck2 * TH1 + btc * Tc - dck2 * IFNg

    # 13. IL-2
    du[13] = bth1ck3 * TH1 - dck3 * IL2
end

end
