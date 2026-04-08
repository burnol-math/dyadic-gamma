=======================
burnolmath/dyadic-gamma
=======================

This page serves as a front-end to
https://gitlab.com/burnolmath/dyadic-gamma.  It is currently quite brief.

.. important::

   If you ended up here by accident and are looking for something serious check:

   1. Brent, Richard P.; Johansson, Fredrik: *A bound for the error term
      in the Brent-McMillan algorithm.* Mathematics of Computation, **84** (295),
      2351–2359. (2015) `doi`__

      __ https://doi.org/10.1090/S0025-5718-2015-02931-7

   2. Yee, A. J.: `Euler-Mascheroni Constant`__ (last updated 2019 attow)

      __ https://www.numberworld.org/digits/EulerGamma/

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

It can be shown that the terms of this alternating series decrease in absolute
value (indeed, this is already the case of the sequence $(e_m/(m+1))_{m\geq1}$), so
their partial sums provide adjacent monotonic series.  However, as long as
some alternative formula is not found for the $e_m$'s, as the recurrence
appears to induce at least a quadratic cost, the series is mainly for fun and
ease of implementation to compute by own modest means from a few dozens to a
few thousands of decimals.

The user interface
==================

The numerical implementation is done using Python's `mpmath <https://mpmath.org>`_.
Examples are given below.

In https://arxiv.org/abs/2603.29998v1 I said that ``ell=4`` or ``ell=5`` were
probably best but actually I obtained here for ``1000`` decimal digits the
shortest running time using ``ell=8``.  So it is the value used by default by the
script.

.. caution::

   Do not use ``ell=1`` with the script, it will take ages to terminate
   already with ``N==1``, as the terminating criterion is (roughly, as it is
   binary not decimal) when the term becomes less than ``10**(-N-3)``, in fact
   here it will compare with ``2**(-15)``, and as the series with $\ell=1$
   behaves about as the alternating series with terms $\log n/n$, the script
   will try to compute more than ``32000`` $e_m$ values, and there is the
   quadratic cost I mentioned.


One can use the script from the command line, as shown below, but, as I was
too lazy to use ``argparse``, only the first two arguments (number of decimal
figures and level) can be passed to the underlying ``gamma_ell()``  function.

The outputs have been hard-wrapped to 50 digits per line, for the convenience
of the reader.

.. code-block:: console

   $ python gammaseries.py 50
   ell is 8
   Last used: m = 24, cm = -2.527e-53
    not used: m = 25, cm = 1.849e-55
   Last digits are 933593992(360)
   gamma rounded to 50 significant figures is (as a string):
   0.57721566490153286060651209008240243104215933593992

   $ python gammaseries.py 50 4
   ell is 4
   Last used: m = 55, cm = 3.021e-52
    not used: m = 56, cm = -3.724e-53
   Last digits are 933593992(363)
   gamma rounded to 50 significant figures is (as a string):
   0.57721566490153286060651209008240243104215933593992

   $ python gammaseries.py 200
   ell is 8
   Last used: m = 95, cm = 7.198e-204
    not used: m = 96, cm = -5.538e-206
   Last digits are 824709491(601)
   gamma rounded to 200 significant figures is (as a string):
   0.57721566490153286060651209008240243104215933593992
     35988057672348848677267776646709369470632917467495
     14631447249807082480960504014486542836224173997644
     92353625350033374293733773767394279259525824709492

   $ python
   < Python 3.13 banner >
   >>> from gammaseries import *
   >>> gamma_ell(400, trunc=True)
   ell is 8
   Last used: m = 189, cm = 2.351e-402
    not used: m = 190, cm = -1.824e-404
   Last digits are 608893312(676)
   gamma truncated to 400 significant figures is (as a string):
   '0.57721566490153286060651209008240243104215933593992
      35988057672348848677267776646709369470632917467495
      14631447249807082480960504014486542836224173997644
      92353625350033374293733773767394279259525824709491
      60087352039481656708532331517766115286211995015079
      84793745085705740029921354786146694029604325421519
      05877553526733139925401296742051375413954911168510
      28079842348775872050384310939973613725530608893312'

References
==========

The first listed paper gives the proof for $\ell\geq2$.  The case $\ell = 1$
is not interesting numerically and requires a separate discussion which was
omitted from the paper.

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

