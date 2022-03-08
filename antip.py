from array import array
from ROOT import TGraphErrors, TFile

output = TFile("antip.root", "recreate")
# http://dx.doi.org/10.1103/PhysRevLett.52.731
# (13-AL-27(AP,NON),,SIG)
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB

pAL = array('f',[474.,548.,626.,715.,782.,881.])
sAL = array('f',[816.,758.,727.,742.,720.,661.])
xAL = array('f',[0. for _ in pAL])
eAL = array('f',[30.,30.,30.,30.,30.,30.])
gAL = TGraphErrors(len(pAL), pAL, sAL, xAL, eAL)
# gAL.Write("13-AL-27")

# http://dx.doi.org/10.1103/PhysRevLett.52.731
# 29-CU-63
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pCU = array('f',[472.,547.,626.,714.,782.,881.])
sCU = array('f',[1220.,1268.,1197.,1217.,1198.,1118.])
xCU = array('f',[0. for _ in pCU])
eCU = array('f',[40.,40.,40.,40.,40.,40.])
gCU = TGraphErrors(len(pCU), pCU, sCU, xCU, eCU)
gCU.Write("29-CU-63")

# http://dx.doi.org/10.1103/PhysRevLett.52.731
# 6-C-12
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pC = array('f',[466.,542.,622.,711.,779.,879.])
sC = array('f',[508.,495.,483.,489.,444.,436.])
xC = array('f',[0. for _ in pC])
eC = array('f',[20.,20.,20.,20.,20.,20.])
gC = TGraphErrors(len(pC), pC, sC, xC, eC)
# gC.Write("6-C-12")

# https://doi.org/10.17182/hepdata.39981
# 4-BE-9
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pBe = array('f',[700,950,1260,1530,1760,2500])
sBe = array('f',[332.29,350.73,325.90,325.79,302.57,282.36])
xBe = array('f',[0. for _ in pBe])
eBe = array('f',[34.70,12.15,8.11,6.52,8.91,5.44])
gBe = TGraphErrors(len(pBe), pBe, sBe, xBe, eBe)
gBe.Write("4-BE-9")

# https://doi.org/10.17182/hepdata.39981
# 6-C-12
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pC = array('f',[700,950,1260,1530,1760,2500])
sC = array('f',[443.50,383.30,380.50,355.11,334.14,323.17])
xC = array('f',[0. for _ in pC])
eC = array('f',[32.85,47.18,9.05,8.24,11.45,11.72])
for i in range(6):
    gC.AddPoint(pC[i],sC[i])
    gC.SetPointError(gC.GetN()-1,xC[i],eC[i])
gC.Write("6-C-12")

# https://doi.org/10.17182/hepdata.39981
# 13-AL-27
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pAl = array('f',[700,950,1260,1530,1760,2500])
sAl = array('f',[805.30,651.97,609.25,609.33,564.05,544.96])
xAl = array('f',[0. for _ in pAl])
eAl = array('f',[33.98,24.19,27.95,13.56,21.96,11.95])
for i in range(6):
    gAL.AddPoint(pAl[i],sAl[i])
    gAL.SetPointError(gAL.GetN()-1,xAl[i],eAl[i])
gAL.Write("13-AL-27")

# https://doi.org/10.17182/hepdata.39981
# 26-FE-56
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pFe = array('f',[700,950,1260,1530,1760,2500])
sFe = array('f',[1206.30,1049.00,1010.30,1034.40,986.34,914.65])
xFe = array('f',[0. for _ in pFe])
eFe = array('f',[134.36,59.75,24.04,21.89,33.09,51.05])
gFe = TGraphErrors(len(pFe), pFe, sFe, xFe, eFe)
gFe.Write("26-FE-56")

# https://doi.org/10.17182/hepdata.39981
# 28-CU-64
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pCu = array('f',[700,950,1260,1530,1760,2500])
sCu = array('f',[1302.40,1115.20,1046.10,1135.20,1010.20,966.24])
xCu = array('f',[0. for _ in pCu])
eCu = array('f',[71.47,34.23,30.54,89.72,184.08,19.62])
gCu = TGraphErrors(len(pCu), pCu, sCu, xCu, eCu)
gCu.Write("28-CU-64")

# https://doi.org/10.17182/hepdata.39981
# 48-CD-112
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pCd = array('f',[700,950,1260,1530,1760,2500])
sCd = array('f',[1831.00,1653.70,1504.50,1562.90,1491.40,1441.10])
xCd = array('f',[0. for _ in pCd])
eCd = array('f',[201.57,76.58,41.92,34.71,59.71,30.76])
gCd = TGraphErrors(len(pCd), pCd, sCd, xCd, eCd)
gCd.Write("48-CD-112")

# https://doi.org/10.17182/hepdata.39981
# 82-PB-207
# MOM         DATA        DATA-ERR
# MEV/C       MB          MB
pPb = array('f',[700,950,1260,1530,1760,2500])
sPb = array('f',[3424.70,2685.50,2463.10,2425.80,2165.20,2120.30])
xPb = array('f',[0. for _ in pPb])
ePb = array('f',[425.11,64.04,109.03,54.33,95.48,71.65])
gPb = TGraphErrors(len(pPb), pPb, sPb, xPb, ePb)
gPb.Write("82-PB-207")