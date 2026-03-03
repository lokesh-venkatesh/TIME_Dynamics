import numpy as np

# log_eps = 1e-12
log_eps = 0

def rhs(t, y, params):
    y = np.maximum(y, log_eps)
    S, SR, C, CR, M1, M2, TH1, TH2, TC, Treg, IL10, IFNgamma, IL2 = y
    
    # Unpack params (matching order from file)
    betaM2, betaTc, betaTh1CK2, betaTh1CK3, betaTh2, betaTreg, gammaC, gammaCR, gammaM1, gammaM2, gammaS, gammaTc, gammaTh1, gammaTh2, gammaTreg, deltaC, deltaCk1, deltaCk2, deltaCk3, deltaCR, deltaM1, deltaM2, deltaS, deltaTc, deltaTh1, deltaTh2, deltaTreg, lambdaM1, lambdaM2, lambdaTc1, lambdaTc2, lambdaTc3, lambdaTc4, lambdaTh1, lambdaTh2, lambdaTreg2, muC1, muC2, muS, muSR, muTcS, muTcTreg, muTh1Ck1, muTh1Ck3, muTregCk1, Cmax, CRmax, k1, k11, k2, k3, k4, k5, k6, k8, k9, ktc1, ktc2, ktc3, ktc4, mC, mS, p1, p2, r1, r2, tck, muM1Ck2, muM2Ck1, k7, k10 = params
    
    dS = (gammaS*((1-mS)*(1-p1-p2)))*S - (deltaS+(p2*gammaS)+gammaS*(mS*p1/2))*S - ((muS*S*IFNgamma)/(IFNgamma+k1)) - ((tck*S*TC)/(ktc1+TC))
    dSR = (gammaS*(1-p1-p2) - (deltaS+(p2*gammaS)))*SR + mS*gammaS*(1-p1/2-p2)*S - ((muSR*SR*IFNgamma)/(k2+IFNgamma)) - ((tck*SR*TC)/(ktc2+TC))
    dC = gammaC*(1-mC)*np.log((Cmax+log_eps)/(C+r1+log_eps))*C + gammaS*(p1+p2)*S - deltaC*C - mC*gammaC*C + (muC1*C*IL10)/(IL10+k3) - (muC2*C*IFNgamma)/(IFNgamma+k4) - (tck*C*TC)/(ktc3+TC)
    dCR = gammaC*CR*np.log((CRmax+log_eps)/(CR+r2+log_eps)) + gammaS*SR*(p1+p2) + mC*gammaC*C - deltaCR*CR + (muC1*CR*IL10)/(IL10+k5) - (muC2*CR*IFNgamma)/(IFNgamma+k6) - (tck*CR*TC)/(ktc4+TC)
    dM1 = gammaM1*M1*((C+CR)/(M1+lambdaM1)) - deltaM1*M1 + ((muM1Ck2*M1*IFNgamma)/(IFNgamma+k7))
    dM2 = gammaM2*M2*((C+CR)/(M2+lambdaM2)) - deltaM2*M2 + ((muM2Ck1*M2*IL10)/(IL10+k10))
    dTH1 = gammaTh1*((TH1*M1)/(lambdaTh1+TH1)) - deltaTh1*TH1 - ((muTh1Ck1*IL10*TH1)/(IL10+k8)) + ((muTh1Ck3*IL2*TH1)/(IL2+k9))
    dTH2 = gammaTh2*((TH2*M2)/(lambdaTh2+TH2)) - deltaTh2*TH2
    dTC = gammaTc*TC*((C+CR)/(TC+lambdaTc1)) + gammaTc*((TC*TH1)/(TC+lambdaTc4)) - muTcS*TC*((S+SR)/(TC+lambdaTc2)) - deltaTc*TC - muTcTreg*TC*((Treg)/(lambdaTc3+Treg))
    dTreg = gammaTreg*((Treg*M2)/(Treg+lambdaTreg2)) - deltaTreg*Treg + ((muTregCk1*IL10*Treg)/(Treg+k11))
    dIL10 = betaM2*M2 - deltaCk1*IL10 + betaTreg*Treg + betaTh2*TH2
    dIFNgamma = betaTh1CK2*TH1 + betaTc*TC - deltaCk2*IFNgamma
    dIL2 = betaTh1CK3*TH1 - deltaCk3*IL2
    
    return np.array([dS, dSR, dC, dCR, dM1, dM2, dTH1, dTH2, dTC, dTreg, dIL10, dIFNgamma, dIL2])
