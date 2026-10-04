// Builds predicates.c, unmodified, from
// https://www.cs.cmu.edu/afs/cs/project/quake/public/code/predicates.c
// (SHA-256 f8662c3f407d1c1c5dcd4dd49ea8b8ddd801a71827d65734206ddc3741792029).
//
// Its error analysis assumes every operation is rounded, so contraction into
// fused multiply-adds is off.

#pragma STDC FP_CONTRACT OFF
#pragma clang diagnostic ignored "-Wdeprecated-non-prototype"

#include "shewchuk.h"
#include "predicates.c"
