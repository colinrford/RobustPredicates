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

#endif
