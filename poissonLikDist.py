import numpy as np

#Measuring Poisson likelihood ratio distance between two normalized count vectors
#In contrary to normal Euclidean distances, this measure considers uncertainty
#of measurements so lower counts gives higher distances
#Note: this is not a metric (doesn't satisfy triangle inequality)


def poissonLikDist(a,b):
    asum = np.sum(a)
    bsum = np.sum(b)


    #MLE expectation of a and b assuming they are drawn from the same distribution
    Ea  = (a+b)*asum/(asum+bsum)
    Eb  = Ea*bsum/asum
    #expectation based on Poisson approximation
    if np.any(Ea==0):
        Ea[Ea==0] =  -np.log(np.mean(a==0))
    if np.any(a==0):
        a[a==0] =  -np.log(np.mean(a==0))
    if np.any(Eb==0):
        Eb[Eb==0] =  -np.log(np.mean(b==0))
    if np.any(b==0):
        b[b==0] =  -np.log(np.mean(b==0))

    logp_unnormalized_adist  = -a + Ea + a*(np.log(a)-np.log(Ea))
    logp_unnormalized_bdist  = -b + Eb + b*(np.log(b)-np.log(Eb))
    return np.sum(logp_unnormalized_adist + logp_unnormalized_bdist)

def poissonLikDistmv(a,b):
    #a is matrix, b is a row vector. Note implementation is different from
    #poissonLikDist for handling zeros
    asum = np.sum(a, axis=1)
    bsum = np.sum(b)

    #MLE expectation of a and b assuming they are drawn from the same distribution
    Eb  = (a+b[np.newaxis,:])*(bsum/(asum+bsum))[:,np.newaxis]
    Ea  = Eb*asum[:,np.newaxis]/bsum

    #expectation based on Poisson approximation.
    if np.any(a==0):
        zeroexp =  -np.log(np.mean(a==0, axis=1))
        a[:,np.isfinite(zeroexp)] = a[:,np.isfinite(zeroexp)] + zeroexp[np.isfinite(zeroexp)]
        Ea[:,np.isfinite(zeroexp)] = Ea[:,np.isfinite(zeroexp)] + zeroexp[np.isfinite(zeroexp)]
    if np.any(b==0):
        zeroexp =  -np.log(np.mean(b==0, axis=1))
        b[:,np.isfinite(zeroexp)] = b[:,np.isfinite(zeroexp)] + zeroexp[np.isfinite(zeroexp)]
        Eb[:,np.isfinite(zeroexp)] = Eb[:,np.isfinite(zeroexp)] + zeroexp[np.isfinite(zeroexp)]

    logp_unnormalized_adist  = -a + Ea + a*(np.log(a)-np.log(Ea))
    logp_unnormalized_bdist  = -b + Eb + b*(np.log(b)-np.log(Eb))
    return np.sum(logp_unnormalized_adist + logp_unnormalized_bdist, axis=1)

def test_example():
    print poissonLikDist(np.asarray([1.,10.]),np.asarray([1.,10.]))
    print poissonLikDist(np.asarray([1.,10.]),np.asarray([0.,10.]))
    print poissonLikDist(np.asarray([1.,10.]),np.asarray([1.,9.]))
    print poissonLikDist(np.asarray([1.,10.]),np.asarray([1.,5.]))
    print poissonLikDist(np.asarray([1.,10.]),np.asarray([1.,20.]))
    print poissonLikDist(np.asarray([1.,5.]),np.asarray([1.,20.]))

    print poissonLikDist(np.asarray([10.,100.]),np.asarray([1.,10.]))
    print poissonLikDist(np.asarray([10.,100.]),np.asarray([1.,9.]))
    print poissonLikDist(np.asarray([10.,100.]),np.asarray([1.,5.]))
    print poissonLikDist(np.asarray([10.,100.]),np.asarray([1.,20.]))

    print poissonLikDistmv(np.asarray([[1,10],[10.,100.]]),np.asarray([1.,9.]))
    print poissonLikDistmv(np.asarray([[0,10],[10.,100.]]),np.asarray([1.,9.]))
