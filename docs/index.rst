=======================
burnolmath/dyadic-gamma
=======================

This page serves as a front-end to
https://gitlab.com/burnolmath/dyadic-gamma.  It is currently quite brief.

.. contents::

The recurrence
==============

$e_0 = 0$ and

$$e_{m} = \frac{2^{m+1} + \sum_{j=1}^m \binom{m+1}{j} e_{m-j}}{2^{m+1} -2}.$$

The series
==========

For any $\ell\geq1$:

$$\gamma = \sum_{1\leq n<2^{\ell-1}}n^{-1} - (\ell-1)\log 2 + \sum_{m=1}^\infty
\frac{(-1)^{m-1}e_m}{m+1}\sum_{2^{\ell-1}\leq n <2^{\ell}}n^{-m-1}$$

The numerical implementation
============================

It is done using Python's `mpmath <https://mpmath.org>`_.  The recurrence has
a total cost to reach a given number of terms which appears to be at least
quadratic, so the method is not competitive for high precision but probably
fine for up to a few hundreds decimal digits.

In https://arxiv.org/abs/2603.29998 I said that ``ell=4`` or ``ell=5`` were
probably best but actually I obtained here for ``1000`` decimal digits the
fastest run time using ``ell=8``.  So it is the value used by default by the
script.

The script is not designed for returning a value for further computations as
``mpmath`` naturally has its own more efficient lazy constant ``euler``.  It
simply prints out the computed value and some extra info.

.. code-block:: console

   $ python gammaseries.py 50
   ell is 8
   Last used: m = 24, |cm| = 2.527e-53
    not used: m = 25, |cm| = 1.849e-55
   gamma = 0.57721566490153286060651209008240243104215933593992(360)

   $ python gammaseries.py 50 4
   ell is 4
   Last used: m = 56, |cm| = 3.724e-53
    not used: m = 57, |cm| = 4.592e-54
   gamma = 0.57721566490153286060651209008240243104215933593992(359)

   $ python gammaseries.py 200
   ell is 8
   Last used: m = 95, |cm| = 7.198e-204
    not used: m = 96, |cm| = 5.538e-206
   gamma = 0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491(601)

References
==========

The first listed paper gives the proof for $\ell\geq2$.  The case $\ell = 1$
is not interesting numerically and requires a separate discussion which was
omitted from the paper. (Do not use ``ell=1`` with the script, it will take
ages to terminate already with ``N==2`` (untested)).

The proof relies on the result of the third paper which itself depends on a
result proven in the fourth paper.  The second paper is quoted only to
increase my citation index.

- Some geometric series for Euler's constant,
  https://arxiv.org/abs/2603.29998

- On the analytic continuation of Dirichlet series with missing digits,
  https://arxiv.org/abs/2602.19727

- Some series representing the eta function for $\Re s>0$,
  https://arxiv.org/abs/2602.05511

- Some series representing the zeta function for $\Re s>1$,
  https://arxiv.org/abs/2601.23158

