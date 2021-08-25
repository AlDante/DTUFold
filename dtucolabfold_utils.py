##########################################################################
##########################################################################
import numpy as np


def kabsch(a, b, weights=None, return_v=False):
    """ Calculate optimal rotation matrix using Kabsch algorithm

    Calculate optimal rotation matrix between two paired sets of points using Kabsch algorithm.
    Optimial in this case means minimum RMSD (root mean square deviation).
    https://en.wikipedia.org/wiki/Kabsch_algorithm
    https://www.pymolwiki.org/index.php/Kabsch

    :param a: first set of points
    :param b: second set of points
    :param weights: weights of chemical elements
    :param return_v: True if vector is to be returned, false to return rotation matrix
    :type return_v: bool
    :return:
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if weights is None:
        weights = np.ones(len(b))
    else:
        weights = np.asarray(weights)

    # B is the covariance matrix, calculated using Einstein summation convention
    B = np.einsum('ji,jk->ik', weights[:, None] * a, b)

    # Singular value decomposition. u * vh is the rotation matrix.
    # This beautiful step provides the answer.  u and vh are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    u, s, vh = np.linalg.svd(B)

    # Convert to right-handed coordinates if necessary.
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation. u and vh are orthonormal, so their det's
    # are +/-1.
    if np.linalg.det(u @ vh) < 0: u[:, -1] = -u[:, -1]

    # return either vector or rotation matrix
    if return_v:
        return u
    else:
        return u @ vh
