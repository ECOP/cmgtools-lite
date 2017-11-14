#include <iostream>

// the final bit: 58 bins
int SRcorr(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=  0 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100                           ) return offset + 12;
        if(mT >= 100 && mT < 160                           ) return offset + 13;
        if(mT >= 160                                       ) return offset + 14;
    }
    if(mll >=  75 && mll < 105 && ht >=  0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return 0          ; // this is the WZ CR !
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 15;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 16;
        if(mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 17;
        if(mT >=   0 && mT < 100 && met >= 250             ) return offset + 18;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 19;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 20;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 21;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 22;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 23;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 24;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 25;
        if(mT >= 160 &&             met >= 200             ) return offset + 26;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 27;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 28;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 29;
        if(mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 30;
        if(mT >=   0 && mT < 100 && met >= 250             ) return offset + 31;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 32;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 33;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 34;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 35;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 36;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 37;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 38;
        if(mT >= 160 &&             met >= 200             ) return offset + 39;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 150) return offset + 40;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 250) return offset + 41;
        if(mT >=   0 && mT < 100 && met >= 250 && met < 350) return offset + 42;
        if(mT >=   0 && mT < 100 && met >= 350             ) return offset + 43;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 44;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 45;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 46;
        if(mT >= 100 && mT < 160 && met >= 200 && met < 250) return offset + 47;
        if(mT >= 100 && mT < 160 && met >= 250 && met < 300) return offset + 48;
        if(mT >= 100 && mT < 160 && met >= 300             ) return offset + 49;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 50;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 51;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 52;
        if(mT >= 160 &&             met >= 200 && met < 250) return offset + 53;
        if(mT >= 160 &&             met >= 250 && met < 300) return offset + 54;
        if(mT >= 160 &&             met >= 300             ) return offset + 55;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100                           ) return offset + 56;
        if(mT >= 100 && mT < 160                           ) return offset + 57;
        if(mT >= 160                                       ) return offset + 58;
    } 
    return 0;
}



// optimization stuff below
// =============================================================================
// =============================================================================

int SRwz(float mT2L, float mT2T, float mll, float mT, float met, int offset = 0) {

    if(mll >= 75 && mll < 105) {
        if (mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 15;
        if (mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 16;
        if (mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 17;
        if (mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 18;
        if (mT >=   0 && mT < 100 && met >= 250 && met < 400) return offset + 19;
        if (mT >=   0 && mT < 100 && met >= 400 && met < 550) return offset + 20;
        if (mT >=   0 && mT < 100 && met >= 550             ) return offset + 21;
        if (mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 22;
        if (mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 23;
        if (mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 24;
        if (mT >= 100 && mT < 160 && met >= 200             ) return offset + 25;
        if (mT >= 160 &&             met >=   0 && met < 100) return offset + 26;
        if (mT >= 160 &&             met >= 100 && met < 150) return offset + 27;
        if (mT >= 160 &&             met >= 150 && met < 200) return offset + 28;
        if (mT >= 160 &&             met >= 200 && met < 250) return offset + 29;
        if (mT >= 160 &&             met >= 250 && met < 400) return offset + 30;
        if (mT >= 160 &&             met >= 400             ) return offset + 31;
    }
    return 0;
}


int SRwzMll(float mll, float mT, float met, int offset = 0) {

    if(mll >= 0 && mll < 75) {
        if (mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if (mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if (mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if (mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if (mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if (mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if (mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset +  7;
        if (mT >= 100 && mT < 160 && met >= 200             ) return offset +  8;
        if (mT >= 160 &&             met >=   0 && met < 100) return offset +  9;
        if (mT >= 160 &&             met >= 100 && met < 150) return offset + 10;
        if (mT >= 160 &&             met >= 150 && met < 200) return offset + 11;
        if (mT >= 160 &&             met >= 200             ) return offset + 12;
    }
    if(mll >= 75 && mll < 105) {
        if (mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 13;
        if (mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 14;
        if (mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 15;
        if (mT >=   0 && mT < 100 && met >= 200             ) return offset + 16;
        if (mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 17;
        if (mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 18;
        if (mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 19;
        if (mT >= 100 && mT < 160 && met >= 200             ) return offset + 20;
        if (mT >= 160 &&             met >=   0 && met < 100) return offset + 21;
        if (mT >= 160 &&             met >= 100 && met < 150) return offset + 22;
        if (mT >= 160 &&             met >= 150 && met < 200) return offset + 23;
        if (mT >= 160 &&             met >= 200             ) return offset + 24;
    }
    if(mll >= 105) {
        if (mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 25;
        if (mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 26;
        if (mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 27;
        if (mT >=   0 && mT < 100 && met >= 200             ) return offset + 28;
        if (mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 29;
        if (mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 30;
        if (mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 31;
        if (mT >= 100 && mT < 160 && met >= 200             ) return offset + 32;
        if (mT >= 160 &&             met >=   0 && met < 100) return offset + 33;
        if (mT >= 160 &&             met >= 100 && met < 150) return offset + 34;
        if (mT >= 160 &&             met >= 150 && met < 200) return offset + 35;
        if (mT >= 160 &&             met >= 200             ) return offset + 36;
    }
    return 0;
}

int SRwz2(float lep1pt, float mT2L, float mT2T, float mll, float mT, float met, int offset = 0) {
    if(lep1pt < 200) 
        return SRwzMll(mll, mT, met, 0);
    else
        return SRwzMll(mll, mT, met, 36);
}

int SRwz3fine(float mll, float mT, float met, float ht, int offset = 0) {
    if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
    if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
    if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
    if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
    if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
    if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
    if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset +  7;
    if(mT >= 100 && mT < 160 && met >= 200             ) return offset +  8;
    if(mT >= 160 &&             met >=   0 && met < 100) return offset +  9;
    if(mT >= 160 &&             met >= 100 && met < 150) return offset + 10;
    if(mT >= 160 &&             met >= 150 && met < 200) return offset + 11;
    if(mT >= 160 &&             met >= 200             ) return offset + 12;
    return 0;
}

int SRwz3finer(float mll, float mT, float met, float ht, int offset = 0) {
    if(mT >=   0 && mT <  40 && met >=   0 && met < 100) return offset +  1;
    if(mT >=   0 && mT <  40 && met >= 100 && met < 150) return offset +  2;
    if(mT >=   0 && mT <  40 && met >= 150 && met < 200) return offset +  3;
    if(mT >=   0 && mT <  40 && met >= 200             ) return offset +  4;
    if(mT >=  40 && mT < 100 && met >=   0 && met < 100) return offset +  5;
    if(mT >=  40 && mT < 100 && met >= 100 && met < 150) return offset +  6;
    if(mT >=  40 && mT < 100 && met >= 150 && met < 200) return offset +  7;
    if(mT >=  40 && mT < 100 && met >= 200             ) return offset +  8;
    if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  9;
    if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 10;
    if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 11;
    if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 12;
    if(mT >= 160 &&             met >=   0 && met < 100) return offset + 13;
    if(mT >= 160 &&             met >= 100 && met < 150) return offset + 14;
    if(mT >= 160 &&             met >= 150 && met < 200) return offset + 15;
    if(mT >= 160 &&             met >= 200             ) return offset + 16;
    return 0;
}


int SRwz3(float mll, float mT, float met, float ht, int offset = 0) {

    if(ht >=   0 && ht <  50) return SRwz3finer(mll, mT, met, ht,  0);
    if(ht >=  50 && ht < 150) return SRwz3finer(mll, mT, met, ht, 16);
    if(ht >= 150 && ht < 300) return SRwz3finer(mll, mT, met, ht, 32);
    if(ht >= 300            ) return SRwz3finer(mll, mT, met, ht, 48);
    return 0;
}

int SRwz4(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=   0 && ht < 100) return SRwz3fine(mll, mT, met, ht,  0);
    if(              mll <  75 && ht >= 100 && ht < 200) return SRwz3fine(mll, mT, met, ht, 12);
    if(              mll <  75 && ht >= 200            ) return SRwz3fine(mll, mT, met, ht, 24);
    if(mll >=  75 && mll < 105 && ht >=   0 && ht < 100) return SRwz3fine(mll, mT, met, ht, 36);
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) return SRwz3fine(mll, mT, met, ht, 48);
    if(mll >=  75 && mll < 105 && ht >= 200            ) return SRwz3fine(mll, mT, met, ht, 60);
    if(mll >= 105              && ht >=   0 && ht < 100) return SRwz3fine(mll, mT, met, ht, 72);
    if(mll >= 105              && ht >= 100 && ht < 200) return SRwz3fine(mll, mT, met, ht, 84);
    if(mll >= 105              && ht >= 200            ) return SRwz3fine(mll, mT, met, ht, 96);
    return 0;
}

int SRwz4p1(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=   0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 12;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 13;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 14;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 15;
        if(mT >= 100 &&             met >=   0 && met < 100) return offset + 16;
        if(mT >= 100 &&             met >= 100 && met < 150) return offset + 17;
        if(mT >= 100 &&             met >= 150             ) return offset + 18;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 19;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 20;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 21;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 22;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 23;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 24;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset + 25;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 26;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 27;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 28;
        if(mT >= 160 &&             met >= 200             ) return offset + 29;
    }
    if(mll >=  75 && mll < 105 && ht >=   0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 30;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 31;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 32;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 33;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 34;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 35;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 36;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 37;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 38;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 39;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 40;
        if(mT >= 160 &&             met >= 200             ) return offset + 41;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 42;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 43;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 44;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 45;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 46;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 47;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 48;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 49;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 50;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 51;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 52;
        if(mT >= 160 &&             met >= 200             ) return offset + 53;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 54;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 55;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 56;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 57;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 58;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 59;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 60;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 61;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 62;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 63;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 64;
        if(mT >= 160 &&             met >= 200             ) return offset + 65;
    }
    if(mll >= 105              && ht >=   0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 66;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 67;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 68;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 69;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 70;
        if(mT >= 160 &&             met >= 100             ) return offset + 71;
    } 
    if(mll >= 105              && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 72;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 73;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 74;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 75;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 76;
        if(mT >= 160 &&             met >= 100             ) return offset + 77;
    } 
    if(mll >= 105              && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 78;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 79;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 80;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 81;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 82;
        if(mT >= 160 &&             met >= 100             ) return offset + 83;
    } 
    return 0;
}

int SRwz4p2(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=   0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 12;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 13;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 14;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 15;
        if(mT >= 100 &&             met >=   0 && met < 100) return offset + 16;
        if(mT >= 100 &&             met >= 100 && met < 150) return offset + 17;
        if(mT >= 100 &&             met >= 150             ) return offset + 18;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 19;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 20;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 21;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 22;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 23;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 24;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset + 25;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 26;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 27;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 28;
        if(mT >= 160 &&             met >= 200             ) return offset + 29;
    }
    if(mll >=  75 && mll < 105 && ht >=   0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 30;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 31;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 32;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 33;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 34;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 35;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 36;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 37;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 38;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 39;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 40;
        if(mT >= 160 &&             met >= 200             ) return offset + 41;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 42;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 43;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 44;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 45;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 46;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 47;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 48;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 49;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 50;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 51;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 52;
        if(mT >= 160 &&             met >= 200             ) return offset + 53;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 54;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 55;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 56;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 57;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 58;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 59;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 60;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 61;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 62;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 63;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 64;
        if(mT >= 160 &&             met >= 200             ) return offset + 65;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 66;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 67;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 68;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 69;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 70;
        if(mT >= 160 &&             met >= 100             ) return offset + 71;
    } 
    return 0;
}

int SRwz4p3(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=  30 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 12;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 13;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 14;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 15;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 16;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 17;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset + 18;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 19;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 20;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 21;
        if(mT >= 160 &&             met >= 200             ) return offset + 22;
    }
    if(mll >=  75 && mll < 105 && ht >=  30 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 23;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 24;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 25;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 26;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 27;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 28;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 29;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 30;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 31;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 32;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 33;
        if(mT >= 160 &&             met >= 200             ) return offset + 34;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 35;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 36;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 37;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 38;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 39;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 40;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 41;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 42;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 43;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 44;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 45;
        if(mT >= 160 &&             met >= 200             ) return offset + 46;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 47;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 48;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 49;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 50;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 51;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 52;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 53;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 54;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 55;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 56;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 57;
        if(mT >= 160 &&             met >= 200             ) return offset + 58;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 59;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 60;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 61;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 62;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 63;
        if(mT >= 160 &&             met >= 100             ) return offset + 64;
    } 
    return 0;
}

int SRwz4p4(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=  0 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 12;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 13;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 14;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 15;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 16;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 17;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset + 18;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 19;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 20;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 21;
        if(mT >= 160 &&             met >= 200             ) return offset + 22;
    }
    if(mll >=  75 && mll < 105 && ht >=  0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 23;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 24;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 25;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 26;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 27;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 28;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 29;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 30;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 31;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 32;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 33;
        if(mT >= 160 &&             met >= 200             ) return offset + 34;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 35;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 36;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 37;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 38;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 39;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 40;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 41;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 42;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 43;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 44;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 45;
        if(mT >= 160 &&             met >= 200             ) return offset + 46;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 47;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 48;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 49;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 50;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 51;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 52;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 53;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 54;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 55;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 56;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 57;
        if(mT >= 160 &&             met >= 200             ) return offset + 58;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 59;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 60;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 61;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 62;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 63;
        if(mT >= 160 &&             met >= 100             ) return offset + 64;
    } 
    return 0;
}

int SRwz4p5(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=  0 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100                           ) return offset + 12;
        if(mT >= 100 && mT < 160                           ) return offset + 13;
        if(mT >= 160                                       ) return offset + 14;
    }
    if(mll >=  75 && mll < 105 && ht >=  0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 15;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 16;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 17;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 18;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 19;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 20;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 21;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 22;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 23;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 24;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 25;
        if(mT >= 160 &&             met >= 200             ) return offset + 26;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 27;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 28;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 29;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 30;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 31;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 32;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 33;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 34;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 35;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 36;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 37;
        if(mT >= 160 &&             met >= 200             ) return offset + 38;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 39;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 40;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 41;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset + 42;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 43;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 44;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 45;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 46;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 47;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 48;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 49;
        if(mT >= 160 &&             met >= 200             ) return offset + 50;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 51;
        if(mT >=   0 && mT < 100 && met >= 100             ) return offset + 52;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 53;
        if(mT >= 100 && mT < 160 && met >= 100             ) return offset + 54;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 55;
        if(mT >= 160 &&             met >= 100             ) return offset + 56;
    } 
    return 0;
}

int SRwz4p6(float mll, float mT, float met, float ht, int offset = 0) {

    if(              mll <  75 && ht >=  0 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset +  1;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset +  2;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset +  3;
        if(mT >=   0 && mT < 100 && met >= 200             ) return offset +  4;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset +  5;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset +  6;
        if(mT >= 100 && mT < 160 && met >= 150             ) return offset +  7;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset +  8;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset +  9;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 10;
        if(mT >= 160 &&             met >= 200             ) return offset + 11;
    }
    if(              mll <  75 && ht >= 200            ) {
        if(mT >=   0 && mT < 100                           ) return offset + 12;
        if(mT >= 100 && mT < 160                           ) return offset + 13;
        if(mT >= 160                                       ) return offset + 14;
    }
    if(mll >=  75 && mll < 105 && ht >=  0 && ht < 100) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 15;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 16;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 17;
        if(mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 18;
        if(mT >=   0 && mT < 100 && met >= 300             ) return offset + 19;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 20;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 21;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 22;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 23;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 24;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 25;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 26;
        if(mT >= 160 &&             met >= 200             ) return offset + 27;
    }
    if(mll >=  75 && mll < 105 && ht >= 100 && ht < 200) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 28;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 29;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 30;
        if(mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 31;
        if(mT >=   0 && mT < 100 && met >= 250             ) return offset + 32;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 33;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 34;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 35;
        if(mT >= 100 && mT < 160 && met >= 200             ) return offset + 36;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 37;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 38;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 39;
        if(mT >= 160 &&             met >= 200             ) return offset + 40;
    }
    if(mll >=  75 && mll < 105 && ht >= 200            ) {
        if(mT >=   0 && mT < 100 && met >=   0 && met < 100) return offset + 41;
        if(mT >=   0 && mT < 100 && met >= 100 && met < 150) return offset + 42;
        if(mT >=   0 && mT < 100 && met >= 150 && met < 200) return offset + 43;
        if(mT >=   0 && mT < 100 && met >= 200 && met < 250) return offset + 44;
        if(mT >=   0 && mT < 100 && met >= 250 && met < 300) return offset + 45;
        if(mT >=   0 && mT < 100 && met >= 300 && met < 350) return offset + 46;
        if(mT >=   0 && mT < 100 && met >= 350             ) return offset + 47;
        if(mT >= 100 && mT < 160 && met >=   0 && met < 100) return offset + 48;
        if(mT >= 100 && mT < 160 && met >= 100 && met < 150) return offset + 49;
        if(mT >= 100 && mT < 160 && met >= 150 && met < 200) return offset + 50;
        if(mT >= 100 && mT < 160 && met >= 200 && met < 250) return offset + 51;
        if(mT >= 100 && mT < 160 && met >= 250 && met < 300) return offset + 52;
        if(mT >= 100 && mT < 160 && met >= 300             ) return offset + 53;
        if(mT >= 160 &&             met >=   0 && met < 100) return offset + 54;
        if(mT >= 160 &&             met >= 100 && met < 150) return offset + 55;
        if(mT >= 160 &&             met >= 150 && met < 200) return offset + 56;
        if(mT >= 160 &&             met >= 200 && met < 250) return offset + 57;
        if(mT >= 160 &&             met >= 250 && met < 300) return offset + 58;
        if(mT >= 160 &&             met >= 300             ) return offset + 59;
    }
    if(mll >= 105                                      ) {
        if(mT >=   0 && mT < 100                           ) return offset + 60;
        if(mT >= 100 && mT < 160                           ) return offset + 61;
        if(mT >= 160                                       ) return offset + 62;
    } 
    return 0;
}


// SSR
// =============================================================================
// =============================================================================

int SuperSigcorr1(float mll, float mT, float met, float ht, int offset = 0) {
    if(mT>100 && met>200 && ht>200) return 1;
    return 0;
}
int SuperSigcorr2(float mll, float mT, float met, float ht, int offset = 0) {
    if(mT>100 && met>200) return 1;
    return 0;
}
int SuperSigcorr3(float mll, float mT, float met, float ht, int offset = 0) {
    if(met>200) return 1;
    return 0;
}
int SuperSigcorr4(float mll, float mT, float met, float ht, int offset = 0) {
    if(mll>86 && mll<96 && met>200) return 1;
    return 0;
}
int SuperSigcorr5(float mll, float mT, float met, float ht, int offset = 0) {
    if(mll>86 && mll<96 && met>300 && ht>200) return 1;
    return 0;
}
int SuperSigcorr6(float mll, float mT, float met, float ht, int offset = 0) {
    if(mll>86 && mll<96 && met>300) return 1;
    return 0;
}



void functionsEWKcorr() {}
