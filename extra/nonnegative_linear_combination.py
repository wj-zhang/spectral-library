import numpy as np
from cvxopt import solvers, matrix


def nlc(C, d, options={'show_progress': False}):
    C = numpy_to_cvxopt_matrix(C)
    d = numpy_to_cvxopt_matrix(d)
    Q = C.T * C
    # try:
    q = - d.T * C
    # except TypeError:
    #     a = 1

    n = C.size[1]
    A = -np.eye(n)
    b = np.zeros(n)
    A = numpy_to_cvxopt_matrix(A)
    b = numpy_to_cvxopt_matrix(b)

    if options:
        for k, v in options.items():
            solvers.options[k] = v

    ret = solvers.qp(Q, q.T, A, b, None, None, None, None)
    w = np.array(ret['x'].T).squeeze()

    if w.shape == ():
        w = np.array([w])
    return w


def numpy_to_cvxopt_matrix(A):
    if A.ndim == 1:
        return matrix(A, (A.shape[0], 1), 'd')
    else:
        return matrix(A, A.shape, 'd')


def optimize_weights_given_alignment(query, W):
    s = query['intensity array']

    len_eles, len_cands = W.shape
    unligned_its = W[-1, :]
    Wa = np.vstack((W[:-1, :], np.zeros((len_cands, len_cands))))
    Wa[np.arange(len_cands) + len_eles - 1, np.arange(len_cands)] = unligned_its

    P = numpy_to_cvxopt_matrix(np.dot(Wa.T, Wa))
    A = numpy_to_cvxopt_matrix(np.dot(s, Wa[:len_eles-1, :])).T
    b = numpy_to_cvxopt_matrix(np.array([1]))

    n = Wa.shape[1]
    G = numpy_to_cvxopt_matrix(-np.eye(n))
    h = numpy_to_cvxopt_matrix(np.zeros(n))
    q = numpy_to_cvxopt_matrix(np.zeros(n))

    solvers.options['show_progress'] = False
    ret = solvers.qp(P, q, G=G, h=h, A=A, b=b)
    w = np.array(ret['x'].T).squeeze()
    if w.size == 1:
        w = np.array([w])
    return w, Wa
