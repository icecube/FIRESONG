# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays


# Random Number Generation
The way we generate random number according to the observed rate function is like this
1. generate a histogram using observed rate function

2. use py.cumsum to sum the histogram

3. define the cumsum as a list called cfd, and normalize cfd to 1, so cfd is a list of number between 0 to 1

4. This is the tricky part, we generate a list of random number, evenly distributed, and try to add them into 
the "cfd" list at the sorted position, i.e. we put a new entry 1.5 between 1 and 2. Then we obtain the index 
of the new entry using np.searchsorted function. Then, we will get the value of the bin using the index.

5. For a larger value of the observed rate, the cumsum increase faster, so more bin value will be pulled from 
that region. The pulled values of the bin mid-point will be the list of the random number. The range of the 
random number will be the range of the x-variable we put into the observed rate function, which is 0 to 1 in 
this case.

6. This new list of random number can be used as log(1+z), which also goes from about 0 to 1.

The number of bin must be large enough in the first step, to ensure the random number generated is fine 
enough.
