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
gAL.Write("13-AL-27")

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
gC.Write("6-C-12")