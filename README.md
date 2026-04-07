# dyadic series for computing gamma

[[_TOC_]]

## Description

Use either within some interactive session:

```console
$ ipython -i gammaseries.py
<Python/IPython banners>

gamma_ell(N[,ell]) gives gamma to N decimal digits.  Default: ell=8

In [1]: help(gamma_ell)


In [2]: gamma_ell(20)
ell is 8
Last used: m = 10, |cm| = 3.218e-23
 not used: m = 11, |cm| = 2.167e-25
Last digits are 153286060(651)
gamma rounded to 20 significant figures is (as a string):
Out[2]: '0.57721566490153286061'

In [3]: gamma_ell(50)
ell is 8
Last used: m = 24, |cm| = 2.527e-53
 not used: m = 25, |cm| = 1.849e-55
Last digits are 933593992(360)
gamma rounded to 50 significant figures is (as a string):
Out[3]: '0.57721566490153286060651209008240243104215933593992'

In [4]: gamma_ell(100)
ell is 8
Last used: m = 47, |cm| = 2.955e-102
 not used: m = 48, |cm| = 2.233e-104
Last digits are 917467495(147)
gamma rounded to 100 significant figures is (as a string):
Out[4]: '0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495'

In [5]: gamma_ell(1000)
ell is 8
Last used: m = 474, |cm| = 2.412e-1003
 not used: m = 475, |cm| = 1.880e-1005
Last digits are 574739302(373)
gamma rounded to 1000 significant figures is (as a string):
Out[5]: '0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495146314472498070824809605040144865428362241739976449235362535003337429373377376739427925952582470949160087352039481656708532331517766115286211995015079847937450857057400299213547861466940296043254215190587755352673313992540129674205137541395491116851028079842348775872050384310939973613725530608893312676001724795378367592713515772261027349291394079843010341777177808815495706610750101619166334015227893586796549725203621287922655595366962817638879272680132431010476505963703947394957638906572967929601009015125195950922243501409349871228247949747195646976318506676129063811051824197444867836380861749455169892792301877391072945781554316005002182844096053772434203285478367015177394398700302370339518328690001558193988042707411542227819716523011073565833967348717650491941812300040654693142999297779569303100503086303418569803231083691640025892970890985486825777364288253954925873629596133298574739302'
```

or via one-shot calls:

```console
$ python gammaseries.py 50
ell is 8
Last used: m = 24, |cm| = 2.527e-53
 not used: m = 25, |cm| = 1.849e-55
Last digits are 933593992(360)
gamma rounded to 50 significant figures is (as a string):
0.57721566490153286060651209008240243104215933593992
```

The output gives `N` decimal digits and then 3 more within parentheses, as the
examples above show the first omitted term should affect only beyond that and
the parenthesized digits, but carries can of course modify them, as well as
other rounding errors (and we do not even bother summing from smallest to
largest, having kept about 3 more decimal guard digits with the hope it is
enough).

The file [`gamma_10000+3`](gamma_10000+3) contains the value of gamma which
was obtained (after enough time for a long coffee break) from:

```console
$ python gammaseries.py 10000
ell is 8
Last used: m = 4745, |cm| = 4.068e-10004
 not used: m = 4746, |cm| = 3.178e-10006
Last digits are 679858165(552)
gamma rounded to 10000 significant figures is (as a string):
0.577215664901532860606512090082402431042159335939923598805...
```

The algorithm is not (at all) suited to as many decimal digits, due to the
cost of the recurrence (see https://burnolmath.gitlab.io/dyadic-gamma for the
underlying mathematical formula).

## References

- Some geometric series for Euler's constant,
  https://arxiv.org/abs/2603.29998

- On the analytic continuation of Dirichlet series with missing digits,
  https://arxiv.org/abs/2602.19727

- Some series representing the eta function for $\Re s>0$,
  https://arxiv.org/abs/2602.05511

- Some series representing the zeta function for $\Re s>1$,
  https://arxiv.org/abs/2601.23158



## License

The files in this repository are distributed under the
CC-BY-SA 4.0 License.  See [LICENSE](LiCENSE).
