"""CRank Bayesian rank aggregation (Zitnik, Sosic & Leskovec, Nat. Commun. 9, 2544, 2018).

A larger value ranks higher; NaN ranks last. nan_policy="order_free" puts every NaN in the worst bin,
weights a tie block that straddles the positive-class boundary fractionally and breaks exact output
ties by mean per-channel rank.
"""
import numpy as np
from scipy import stats

NAN_POLICIES = ('nan_last', 'order_free')

def _rankdata_nan_last(a, nan_policy='nan_last'):
    if nan_policy not in NAN_POLICIES:
        raise ValueError(f'nan_policy {nan_policy!r} not in {NAN_POLICIES}')
    a = np.asarray(a, dtype=float)
    nan_mask = np.isnan(a)
    if not nan_mask.any():
        return stats.rankdata(a)
    n = a.size
    ranks = np.empty(n, dtype=float)
    finite_idx = np.flatnonzero(~nan_mask)
    nan_idx = np.flatnonzero(nan_mask)
    ranks[finite_idx] = stats.rankdata(a[finite_idx])
    ranks[nan_idx] = float(n) if nan_policy == 'order_free' else np.arange(finite_idx.size + 1, n + 1, dtype=float)
    return ranks

def _smooth(x):
    smooth = x.copy()
    n_changes = 1
    while n_changes != 0:
        prev = smooth.copy()
        for i in range(1, len(x) - 1):
            smooth[i] = np.median(prev[i - 1:i + 2])
        smooth[0] = np.median([prev[0], smooth[1], 3 * smooth[1] - 2 * smooth[2]])
        smooth[-1] = np.median([prev[-1], smooth[-2], 3 * smooth[-2] - 2 * smooth[-3]])
        n_changes = np.sum(smooth != prev)
    return smooth

def _cummax(x):
    y = np.array([np.max(x[:i]) for i in range(1, len(x) + 1)])
    return y

def crank_rra(data, prior=0.2, num_bin=4, num_iter=10, return_all=True, corr_stop=1, nan_policy='nan_last', verbose=False):
    data = np.asarray(data, dtype=float)
    nr, nc = data.shape
    nrp = int(np.floor(nr * prior))
    if verbose:
        print('Nrp: %d' % nrp)
    rank_data = np.zeros(data.shape)
    for i in range(nc):
        rank_data[:, i] = _rankdata_nan_last(-data[:, i], nan_policy) / float(nr)
    bayes_factors = np.zeros((num_bin, nc))
    binned_data = np.ceil(rank_data * num_bin).astype('int')
    bayes_data = np.zeros((nr, nc))
    guess = np.mean(rank_data, 1)
    cprev = 0
    for iter in range(num_iter):
        if corr_stop - cprev < 1e-15:
            if verbose:
                print('Converged!')
            break
        guess_last = guess.copy()
        if nan_policy == 'order_free':
            cut = np.partition(guess, nrp - 1)[nrp - 1]
            better, tied = (guess < cut, guess == cut)
            pos = better.astype(float)
            pos[tied] = (nrp - better.sum()) / tied.sum()
            guess = pos
        else:
            oo = np.argsort(guess)
            guess[oo[:nrp]] = 1.0
            guess[oo[nrp:]] = 0.0
        for i in range(nc):
            for bin in range(1, num_bin + 1):
                tpr = np.sum(guess[binned_data[:, i] <= bin])
                fpr = np.sum((1.0 - guess)[binned_data[:, i] <= bin])
                bayes_factors[bin - 1, i] = np.log((tpr + 1.0) / (fpr + 1.0) / (prior / (1.0 - prior)))
        for i in range(nc):
            bayes_factors[:, i] = _smooth(bayes_factors[:, i])
            bayes_factors[:, i] = _cummax(bayes_factors[:, i][::-1])[::-1]
        for i in range(nc):
            bayes_data[:, i] = bayes_factors[binned_data[:, i] - 1, i]
        guess = stats.rankdata(-np.sum(bayes_data, 1))
        cprev = stats.pearsonr(guess, guess_last)[0]
        if verbose:
            print('Correlation with previous iteration: %f' % cprev)
    if nan_policy == 'order_free':
        s_rank = stats.rankdata(-np.sum(bayes_data, 1), method='dense').astype(np.int64)
        m_rank = stats.rankdata(np.mean(rank_data, 1), method='dense').astype(np.int64)
        guess = stats.rankdata(s_rank * (m_rank.max() + 1) + m_rank)
    if return_all:
        return (guess, bayes_data, bayes_factors)
    else:
        return guess
