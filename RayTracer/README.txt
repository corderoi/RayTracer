Ian Cordero | UID 704058197
07 December 2015
CS 174A | Prof. Diana Ford | UCLA Fall 2015
Dis 1B | TA: Donyang

--

Test Results
------------

Test                                | Result
1.  testBehind.ppm: 600 x 600:        PASS
2.  testShadow.ppm: 600 x 600:        PASS
3.  testDiffuse.ppm: 600 x 600:       PASS
4.  testSpecular.ppm: 600 x 600:      PASS
5.  testAmbient.ppm: 600 x 600:       PASS
6.  testParsing.ppm: 600 x 600:       PASS
7.  testSample.ppm: 600 x 600:        PASS
8.  testReflection.ppm: 600 x 600:    PASS
9.  testIntersection.ppm: 600 x 600:  PASS
10. testImgPlane.ppm: 600 x 600:      PASS
11. testBackground.ppm: 600 x 600:    PASS
12. testIllum.ppm: 600 x 600:         PASS

This raytracer program passed 12/12 tests. Tests are based on close
visual inspection of the results.

Notes: an exact match (compared using the 'diff' program) to the test
images does not seem possible, hence testing using visual inspection
was deemed necessary.

The program can be compiled via any OS X or linux command line with g++:
'g++ raytrace.cpp -o raytrace'.

Further explanation of this program is detailed in the comments in
'raytrace.cpp'.
