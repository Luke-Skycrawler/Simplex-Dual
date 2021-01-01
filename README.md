# Simplex-Dual
an naive, test-rich implementation of simplex algorithm 

unit tests
[ ] 0 solution
[ ] inf
[ ] ranks

---

#### Standard Form for a LP:

$$
\begin{aligned}
\max\quad &z=c_1x_1+c_2x_2+\cdots+c_nx_n\\
{\rm s.t.}\quad& a_{11}x_1+a_{12}x_2+\cdots+a_{1n'}x_{n'}= b_1\ge0\\
& a_{21}x_1+a_{22}x_2+\cdots+a_{2n'}x_{n'}= b_2\ge0\\
&\cdots\\
& a_{m1}x_1+a_{m2}x_2+\cdots+a_{mn'}x_{n'}= b_m\ge0\\
& x_1\ge0\\
& x_2\ge0\\
& \cdots\\
& x_{n'}\ge0\\
\end{aligned}
$$

To transform a input like:

```
n m
c1 c2 … cn
a11 a12 … a1n b1 d1
a21 a22 … a2n b2 d2
…
am1 am2 … amn bm dm
e1 e2 … en
```

into the standard form,

* If `di` = -1, add a nonnegative auxiliary variable to the left.
* If `di` = 1, subtract a nonnegative auxiliary variable to the left.
* If `ej` = -1, flip over `aij` (i = 1, 2, ..., m) to `-aij` so we can keep `xj` nonnegative during the program; but remember to output `-xj`!
* If `ej` = 0, transform `xj` into the difference of 2 nonnegative auxiliary variables `xm-xn` (so 1 more variable only); remember to output `xm-xn` as `xj`!
* If `bi` < 0, flip over both sides of the equation (do this after all things above are done to make code easier)

#### Primal Dual -- Things to Do with a Standard LP:

1. Choose the initial base variables with **some rule**, meanwhile transfer every row conforming to the base variables
2. Enter the loop:
   1. Choose a non-base variable `xj` with `cj` > 0 to become the new base variable with **some rule**, replacing the old base variable `xi` with min(`bi`/`aij`) (with `aij` > 0; if no such `i`, the optimal value is **infinite**)
   2. Update every other row in the simplex table; then check if all `cj` < 0
      * True: done
      * False: continue looping

