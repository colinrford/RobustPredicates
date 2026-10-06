// Shewchuk's predicates.c, for comparison in tests.

#ifndef SHEWCHUK_H
#define SHEWCHUK_H

void exactinit(void) __attribute__((swift_name("exactinitC()")));

double orient2d(double *pa, double *pb, double *pc)
  __attribute__((swift_name("orient2dC(_:_:_:)")));

double orient3d(double *pa, double *pb, double *pc, double *pd)
  __attribute__((swift_name("orient3dC(_:_:_:_:)")));

double incircle(double *pa, double *pb, double *pc, double *pd)
  __attribute__((swift_name("inCircleC(_:_:_:_:)")));

double insphere(double *pa, double *pb, double *pc, double *pd, double *pe)
  __attribute__((swift_name("inSphereC(_:_:_:_:_:)")));

// Stages B to D of orient2d, called when stage A can't decide.

double orient2dadapt(double *pa, double *pb, double *pc, double detsum)
  __attribute__((swift_name("orient2dAdaptC(_:_:_:detsum:)")));

// Stages B to D of incircle.

double incircleadapt(double *pa, double *pb, double *pc, double *pd, double permanent)
  __attribute__((swift_name("inCircleAdaptC(_:_:_:_:permanent:)")));

// Exact without filtering or adaptivity, for timing.

double orient2dslow(double *pa, double *pb, double *pc)
  __attribute__((swift_name("orient2dSlowC(_:_:_:)")));

double orient3dslow(double *pa, double *pb, double *pc, double *pd)
  __attribute__((swift_name("orient3dSlowC(_:_:_:_:)")));

double incircleslow(double *pa, double *pb, double *pc, double *pd)
  __attribute__((swift_name("inCircleSlowC(_:_:_:_:)")));

#endif
