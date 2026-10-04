Test data from A. Nanevski, G. Blelloch and R. Harper, "Automatic Generation of
Staged Geometric Predicates", Higher-Order and Symbolic Computation
16(4):379–400, 2003 (CMU-CS-01-141, 2001), retrieved from
https://www.cs.cmu.edu/afs/cs/project/pscico/pscico/src/arithmetic/compiler1/test/

| File           | Original      |
| -------------- | ------------- |
| `orient2d.txt` | `orient.2d`   |
| `incircle.txt` | `insphere.2d` |
| `orient3d.txt` | `orient.3d`   |
| `insphere.txt` | `insphere.3d` |

Line endings are converted from CRLF to LF; nothing else is changed. Each line
is an index, the coordinates of the points in order, and the sign of the
determinant. Coordinates are uniformly random over Doubles with exponents
between −63 and 63 (§5).
