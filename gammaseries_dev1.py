"""
Computation of gamma for a given target decimal precision and level

Author: Jean-François Burnol
Created: April 2, 2026 v0.1.0
This version: April 18, 2026 v0.3.0dev

change log:

0.3.0dev: BUT IS NOT EFFICIENT (see below)
- efficiency gain: when loop parameter m becomes larger than the
  needed running target binary precision, stop using exact big int
  for binomial coefficients, switch to using mpf type for them as
  well.  Also, avoid big pure powers of two, due to costly
  conversion times to mpf.  (about this, cf
  https://github.com/mpmath/mpmath/issues/1073)

  BUT IN FACT NO GAIN (to the contrary, even):

    The transition to using mpf for binomial coefficients is done
    only after the bulk of the computation has been done,
    especially with the default ell=8 (which seems best choice in
    a wide range around N=1000): happens at about m=Pmax/ell, of a
    total of about Pmax/(ell-1) computed terms.

0.2.1:
- indicate explicitly the sign of the last used
  and first non used coefficients.

0.2.0:
- add silent and trunc booolean optional parameters,
- return a string in place of printing it.

0.1.0: initial version.

TODO: should I return an object of type mpf rather?

© 2026 Jean-François Burnol

This script is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.
Full license text: https://creativecommons.org/licenses/by-sa/4.0/

ANY RE-USE OR PLAGIARIZING BY AN ARTIFICIAL INTELLIGENCE WITHOUT
PROPER ATTRIBUTION IS STRICTLY FORBIDDEN AND WILL GET PUNISHED

"""

__version__ = "0.3.0"
__date__    = "2026-04-18"
__author__  = "Jean-François Burnol"

from math import log, ceil
from mpmath import mp, mpf, nstr

def gamma_ell(N:int,
              ell:int = 8,
              silent:bool = False,
              trunc:bool = False) -> str:
    """Computes gamma with N decimal digits, using level ell

    :param N int: requested number of decimal places
    :param ell int: (optional, defaults to 8) the level ell
    :param silent bool: (optional, defaults to False) whether to print
                        some extra info
    :param trunc bool: (optional, defaults to False) whether to round or
                       truncate the return value
    :rtype: str
    :return: String with N significant figures of the decimal expansion
             of Euler's constant, either rounded or truncated as per trunc
             optional parameter.

    """
    P_target = ceil((N + 1) * log(10, 2))
    P_eps = P_target + min(ell, 10)
    P_max = P_eps + 12

    mp.prec = P_max
    current_prec = P_max
    epsilon = mpf(2) ** -P_eps

    ell_m_one = ell - 1
    two_to_the_ell_m_one = 2 ** ell_m_one
    two_to_the_ell = 2 * two_to_the_ell_m_one

    S = sum(mpf(1)/n for n in range(two_to_the_ell_m_one - 1, 0, -1))
    S -= (ell-1) * mp.log(mpf(2))

    invpowers = [ 1/mpf(n) for n in range(two_to_the_ell - 1,
                                          two_to_the_ell_m_one - 1,
                                          -1)
                 ]

    cms = [ 0 ]
    ems = [ 0 ]

    m = 0
    two_to_the_m_plus_one = 2
    cmsign = -1

    # 0.3.0 reorganizes into two while loop in place of one
    #       (to avoid repeting "if" tests in the context
    #        of a threshold from big int to mpf types)
    done = False
    while m < current_prec:
        m += 1
        two_to_the_m_plus_one *= 2
        cmsign = -cmsign

        mp.prec = current_prec

        em = mpf(0)
        binomj = 1  # binomial coefficients of type int
        for j in range(1, m):  # e0 = 0 de toute façon
            binomj = binomj * (m - j + 2) // j
            em += mpf(binomj) * ems[-j]

        em += two_to_the_m_plus_one
        em /= two_to_the_m_plus_one - 2
        ems.append(em)

        # keep inverse powers from smallest to largest
        i = 0
        for n in range(two_to_the_ell - 1,
                       two_to_the_ell_m_one - 1,
                       -1):
            invpowers[i] /= n
            i += 1

        # would it be more efficient to use invpowers[::-1]
        # and keep in ordering induced from natural numbers?
        cm = em / ( m + 1 ) * sum(x for x in invpowers)

        if cm < epsilon:
            done = True
            break

        cm *= cmsign
        cms.append(cm)

        mp.prec = P_max
        S += cm

        # 30 binary digits looks a bit excessive.
        current_prec = max(30, current_prec - ell_m_one)
        mp.prec = current_prec

    if not done:
        # Strange looking code, but https://github.com/mpmath/mpmath/issues/1073
        # Anyhow done only once so probably does not matter.
        # The point is that from now on we have an mpf type.
        # But the more significant change is for binomial coefficeints, see below.
        two_to_the_m_plus_one = mpf(two_to_the_m_plus_one + 1) - mpf(1)
        if not silent:
            print(f"TRANSITION à {m}. P_max est {P_max}. current_prec est {current_prec}.\n"
                  f"P_max/ell est {P_max/ell:.2f}. P_max/(ell-1) est {P_max/(ell-1):.2f}")
        # m is about P_max/ell.  There are about P_max/(ell-1) steps.
        # The transition happens when the bulk has been done already.
        # Not advantageous.
        while True:
            m += 1
            two_to_the_m_plus_one *= 2
            cmsign = -cmsign

            mp.prec = current_prec

            em = mpf(0)
            binomj = mpf(1)  # m is at least comparable to current precision
                             # and we stop using exact integers.
            for j in range(1, m):
                # should I worry about conversion efficiency of these int's to mpf ?
                binomj = binomj * (m - j + 2) / j
                em += binomj * ems[-j]

            # Should I use some extra guards here, now that I ignore a "-2"
            # in the denominator? This "-2" would not change the denominator
            # anyway if not immediately soon and we have guard digits already.
            em /= two_to_the_m_plus_one
            # So now we add 1 after (MEMO em is about log2(m)) and not 2**(m+1)
            # before division by 2**(m+1)-2
            # Memo: em is of size about log2(m).  The 1 here is NOT AT ALL negligible.
            em += 1
            ems.append(em)

            # keep inverse powers from smallest to largest
            i = 0
            for n in range(two_to_the_ell - 1,
                           two_to_the_ell_m_one - 1,
                           -1):
                invpowers[i] /= n
                i += 1

            # some inefficiency will happen when m+1 is a big power
            # of 2 (see github issue above) but well
            cm = em / ( m + 1 ) * sum(x for x in invpowers)

            if cm < epsilon:
                break

            cm *= cmsign
            cms.append(cm)

            mp.prec = P_max
            S += cm

            # 30 binary digits looks a bit excessive.
            current_prec = max(30, current_prec - ell_m_one)
            mp.prec = current_prec

    mp.prec = P_max

    if not silent:
        print(f"ell is {ell}")
        print(f"Last used: m = {m-1}, cm = {nstr(cms[-1], 4, strip_zeros=False)}")
        print(f" not used: m = {m}, cm = {nstr(cmsign * cm, 4, strip_zeros=False)}")

        # ATTENTION ! nstr supprime trailing zeros par défaut !!
        S_string = nstr(S, N + 3, strip_zeros=False)
        # print last digits + 3 more
        print(f"Last digits are {S_string[-12:-3]}({S_string[-3:]})")

        print(f"gamma {'truncated' if trunc else 'rounded'} "
              f"to {N} significant figures is (as a string):")

    if trunc:
        # we hope this gives correct truncation...
        return nstr(S, N + 6, strip_zeros=False)[:-6]
    else:
        return nstr(S, N, strip_zeros=False)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        print(gamma_ell(int(sys.argv[1]), int(sys.argv[2])))
    elif len(sys.argv) > 1:
        print(gamma_ell(int(sys.argv[1])))
    else:
        print("\ngamma_ell(N[,ell]) gives gamma to N decimal digits. "
              "Default: ell=8.\nIssue help(gamma_ell) for the other arguments.")
