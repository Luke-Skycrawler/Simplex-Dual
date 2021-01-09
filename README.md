# Simplex-Dual
an naive, test-rich implementation of simplex algorithm 

unit tests
[ ] 0 solution
[ ] inf
[ ] ranks

---

### Primal

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

#### Primal -- Things to Do with a Standard LP:

1. Implement the **2-phase method**, where $m$ variables would be added.
2. Choose the artificial variables as the initial base variables. Actually, you just need to **continue your work** in phase 1. (some cuttings are necessary)
3. Enter the loop:
   1. Choose a non-base variable `xj` with `cj` > 0 to become the new base variable with **some rule**, replacing the old base variable `xi` with min(`bi`/`aij`) (with `aij` > 0; if no such `i`, the optimal value is **infinite**)
   2. Update every other row in the simplex table; then check if all `cj` <= 0
      * True: done
      * False: continue looping

### Dual

#### Standard Form for a LP:

$$
\begin{aligned}
\max\quad &z=c_1x_1+c_2x_2+\cdots+c_nx_n\\
&c_1,c_2,\cdots,c_n\le 0\\
{\rm s.t.}\quad& a_{11}x_1+a_{12}x_2+\cdots+a_{1n'}x_{n'}\le b_1\\
& a_{21}x_1+a_{22}x_2+\cdots+a_{2n'}x_{n'}\le b_2\\
&\cdots\\
& a_{m1}x_1+a_{m2}x_2+\cdots+a_{mn'}x_{n'}\le b_m\\
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

* If `di` = 1, flip the signs of coefficients on both sides.
* If `di` = 0, turn such 'A = b' into 'A <= b' and 'A >= b'.
* If `ej` = -1, flip over `aij` (i = 1, 2, ..., m) to `-aij` so we can keep `xj` nonnegative during the program; but remember to output `-xj`!
* If `ej` = 0, transform `xj` into the difference of 2 nonnegative auxiliary variables `xm-xn` (so 1 more variable only); remember to output `xm-xn` as `xj`!
* After all this, if some **`cj` < 0, don't use dual** :blush:

#### Dual -- Things to Do with a Standard LP:

1. No need of 2-phase. We assure that `cj` < 0 at the beginning. And after standardization, we can easily add a slack variable for each row since we've changed '=' and '>=' into '<=' previously.
2. Choose the the slack variables as initial bases; no need of transforming.
3. Enter the loop:
   1. Choose a base variable `xi` with `bi` < 0 (`bi` should be smallest?) out of base, and select the `xj` with $\min_{\overline a_{ij}<0}\cfrac {\sigma_j}{\overline a_{ij}}$ to be the new base variable.
      * If no `aij < 0`, the dual problem has no finite solution <=> the original problem has no solution.
   2. Update every other row in the simplex table; then check if all `bi` >= 0
      * True: go to 3.
      * False: continue looping