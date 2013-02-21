snp_vis
=======

Benjamin VanderSluis
bvander@cs.umn.edu
February 2013

This is a processing sketch for visualizing my 23&me
SNP data. I intend to one day make it an Xscreensaver 
module, if I can ever figure out how...

It is object oriented, and there are currently two
types of snp objects: spheres and helices. They 
are statically initiated, so you get the same results
each time, but adding new shapes, and new locations,
(even randomly chosen) should be trivial. My current
aim is to get all the pieces working.

**BYOG(enome)**

of the form
snip_id<tab>chromosome_number(## | 'X' | 'Y')<tab>SNP

where SNP is two letters from "ATCG" (or a single letter if
you're a guy and SNP resides on a sex chromosome


