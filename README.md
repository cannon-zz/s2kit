# S2kit: A Lite Version of SpharmonicKit

This is a snapshot of a library for Fourier transforms on the 2-sphere by Kostelec P. and Rockmore D. N. written some time around 2004.  I downloaded the code from their website at the University of Dartmouth.  The original package has since disappeared from the internet, so because I had a need to collaborate with people and use this code I posted it to git.ligo.org in April of 2019.  I believe the initial commit is exactly their original code, without any modifications.

I have accumulated some patches since then.  The directory layout has been reorganized and the build system converted to autotools to make packaging easy.  One or two performance improvements and bugs have been fixed.

The original README file stated

>See
>
> s2kit_fx.pdf

which is now in doc/
