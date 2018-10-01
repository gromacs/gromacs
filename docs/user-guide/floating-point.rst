.. _gmx-floating-point:

Floating point arithmetic
=========================

|Gromacs| spends its life doing arithmetic on real numbers, often summing many
millions of them. These real numbers are encoded on computers in so-called
binary floating-point representation. This representation is somewhat like
scientific exponential notation (but uses binary rather than decimal), and is
necessary for the fastest possible speed for calculations. Unfortunately the
laws of algebra only approximately apply to binary floating-point. In part,
this is because some real numbers that are represented simply and exactly in
decimal (like 1/5=0.2) have no exact representation in binary floating-point,
just as 1/3 cannot be represented in decimal. There are many sources you can
find with a search engine that discuss this issue more exhaustively, such as
`Wikipedia <https://en.wikipedia.org/wiki/Floating-point_arithmetic>`__ and
David Goldberg's 1991 paper *What every computer scientist should know about
floating-point arithmetic* (`article <https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html>`__,
`addendum <https://docs.oracle.com/cd/E37069_01/html/E39019/z400228248508.html>`__).
Bruce Dawson also has a written a number of very valuable blog posts on modern
floating-point programming at his
`Random ASCII site <https://randomascii.wordpress.com/category/floating-point/>`__
that are worth reading.

So, the sum of a large number of binary representations of exact decimal
numbers need not equal the expected algebraic or decimal result. Users observe
this phenomenon in sums of partial charges expressed to two decimal places that
sometimes only approximate the integer total charge to which they contribute
(however a deviation in the first decimal place would always be indicative of a
badly-formed topology).  When |Gromacs| has to represent such floating-point
numbers in output, it sometimes uses a computer form of scientific notation
known as E notation. In such notation, a number like -9.999971e-01 is actually
-0.9999971, which is close enough to -1 for purposes of assessing the total
charge of a system.

It is also not appropriate for |Gromacs| to guess to round things, because such
rounding relies on assumptions about the inputs that need not be true. Instead
the user needs to understand how their tools work.
