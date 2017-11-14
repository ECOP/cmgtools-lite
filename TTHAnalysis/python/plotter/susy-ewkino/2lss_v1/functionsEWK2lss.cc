#include <cmath>
#include <map>

float mt_2_2lss(float pt1, float phi1, float pt2, float phi2) {
    return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

float mtwmin(float ptl1, float phil1, float ptl2, float phil2, float met, float metphi){
  float mtw1 = mt_2_2lss(ptl1, phil1, met, metphi);
  float mtw2 = mt_2_2lss(ptl2, phil2, met, metphi);
  return std::min(mtw1,mtw2);
}

float pt_2_2lss(float pt1, float phi1, float pt2, float phi2) {
    phi2 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2), pt2*std::sin(phi2));
}

int SSR_2lss(int nj, float ptl1, float phil1, float ptl2, float phil2, float met, float metphi, int charge){
  float mtw = mtwmin(ptl1, phil1, ptl2, phil2, met, metphi);
  float ptdil = pt_2_2lss(ptl1,phil1,ptl2,phil2);
  int   cha  = charge;

  if      (nj == 0 && mtw > 100 &&                met > 140) return 1;
  else if (nj == 1 && mtw < 100 && ptdil < 100 && met > 200) return 2;
  else return -1;
}

int SSR_2lss_1(int nj, float ptl1, float phil1, float ptl2, float phil2, float met, float metphi, int charge){
  float mtw1 = mt_2_2lss(ptl1,phil1, met, metphi);
  float mtw2 = mt_2_2lss(ptl2,phil2, met, metphi);
  float mtw  = std::min(mtw1,mtw2);
  float ptdil = pt_2_2lss(ptl1,phil1,ptl2,phil2);
  int   cha  = charge;

  if(nj == 0 && mtw > 100 &&                met > 140) return 1;
  return 0;
}

int SSR_2lss_2(int nj, float ptl1, float phil1, float ptl2, float phil2, float met, float metphi, int charge){
  float mtw1 = mt_2_2lss(ptl1,phil1, met, metphi);
  float mtw2 = mt_2_2lss(ptl2,phil2, met, metphi);
  float mtw  = std::min(mtw1,mtw2);
  float ptdil = pt_2_2lss(ptl1,phil1,ptl2,phil2);
  int   cha  = charge;

  if (nj == 1 && mtw < 100 && ptdil < 100 && met > 200) return 1;
  return 0;
}


int SR_ewk_ss2l(int nj, float ptl1, float phil1, float ptl2, float phil2, float met, float metphi, int charge){
  
  float mtw1 = mt_2_2lss(ptl1,phil1, met, metphi);
  float mtw2 = mt_2_2lss(ptl2,phil2, met, metphi);
  float mtw  = std::min(mtw1,mtw2);
  float ptdil = pt_2_2lss(ptl1,phil1,ptl2,phil2);
  int   cha  = charge;

  if      (nj==0 && ptdil<75 && mtw<100 && met<100)            return 1;  //VR

  else if (nj==0 && ptdil<75 && mtw<100 && met>100 && met<150 && cha == +1) return 2;
  else if (nj==0 && ptdil<75 && mtw<100 && met>100 && met<150 && cha == -1) return 3;


  else if (nj==0 && ptdil<75 && mtw<100 && met>150 && met < 200)            return 4;
  else if (nj==0 && ptdil<75 && mtw<100 && met > 200)            			    return 5;

  else if (nj==0 && ptdil>75 && mtw<100 && met<100)            return 6; //VR


  else if (nj==0 && ptdil>75 && mtw<100 && met>100 && met<150 && cha == +1) return 7;
  else if (nj==0 && ptdil>75 && mtw<100 && met>100 && met<150 && cha == -1) return 8;

  else if (nj==0 && ptdil>75 && mtw<100 && met>150 && met < 200)            return 9;
  else if (nj==0 && ptdil>75 && mtw<100 && met>200)				             return 10;

  else if (nj==0 && mtw>100 && met<100)                        return 11;  //VR


  else if (nj==0 && mtw>100 && met>100 && met<150 && cha == +1)             return 12;
  else if (nj==0 && mtw>100 && met>100 && met<150 && cha == -1)             return 13;

  else if (nj==0 && mtw>100 && met>150 && met<200)                          return 14;
  else if (nj==0 && mtw>100 && met>200)                        			    return 15;

  else if (nj==1 && ptdil<75 && mtw<100 && met<100 )            return 16;  //VR

  else if (nj==1 && ptdil<75 && mtw<100 && met>100 && met<150 && cha == +1) return 17;
  else if (nj==1 && ptdil<75 && mtw<100 && met>100 && met<150 && cha == -1) return 18;

  else if (nj==1 && ptdil<75 && mtw<100 && met>150 && met<200)              return 19;
  else if (nj==1 && ptdil<75 && mtw<100 && met>200)         		          return 20;

  else if (nj==1 && ptdil>75 && mtw<100 && met<100)            return 21; //VR


  else if (nj==1 && ptdil>75 && mtw<100 && met>100 && met<150 && cha == +1) return 22;
  else if (nj==1 && ptdil>75 && mtw<100 && met>100 && met<150 && cha == -1) return 23;

  else if (nj==1 && ptdil>75 && mtw<100 && met>150 && met < 200)            return 24;
  else if (nj==1 && ptdil>75 && mtw<100 && met>200) 					          return 25;

  else if (nj==1 && mtw>100 && met<100)                        return 26;  //VR

  else if (nj==1 && mtw>100 && met>100 && met<150 && cha == +1)             return 27;
  else if (nj==1 && mtw>100 && met>100 && met<150 && cha == -1)             return 28;

  else if (nj==1 && mtw>100 && met>150 && met<200)                          return 29;
  else if (nj==1 && mtw>100 && met>200)                   				       return 30;

  return -99;   

} 

int allTight(int th1, int th2) {
	return (th1&&th2);
}

int sameSign(int ch1, int ch2) {
    return (ch1*ch2>0);
}

int isGoodFake2L(float pt, int isTight) {
    if(pt == 0) return 0;
    if(isTight) return 0;
    return 1;
}



