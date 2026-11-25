"""Conservative frozen-coefficient source map for a C1-smoothed B2 jump.

No transport, porosity, Darcy, permeability or viscosity relation is changed.
Linear charging/transition and the square-law aging portions are integrated
analytically until the smoothing band is reached. Inside that band an adaptive
SSP-RK3 map with step doubling advances the same autonomous scalar ODE.
"""
import numpy as np

DIAGNOSTICS = {}


def reset_diagnostics():
    DIAGNOSTICS.clear()
    DIAGNOSTICS.update(calls=0, smooth_values=0, accepted_substeps=0,
                       rejected_substeps=0, maximum_substeps_per_value=0,
                       maximum_accepted_estimated_error=0.,
                       minimum_internal_deposit=np.inf)


reset_diagnostics()


def smooth_weight(B, T, delta):
    s = (B - (T - delta)) / (2 * delta)
    z = np.minimum(1., np.maximum(0., s))
    return z * z * (3 - 2 * z)


def net_source(B, c, alpha, beta, p):
    lo = np.asarray(p.B1)[:, None]
    T = np.asarray(p.B2)[:, None]
    cap = np.asarray(p.cap)[:, None]
    if p.smoothing_fraction == 0:
        weight = (B > T).astype(float)
    else:
        weight = smooth_weight(B, T, T * p.smoothing_fraction)
    quotient = np.divide(cap, B, out=np.zeros_like(B), where=B > 0)
    q = alpha * c * (1 + weight * (quotient - 1)) - beta * B
    q = np.where(B < lo, np.asarray(p.charging)[:, None] * alpha * c, q)
    q = np.where((B == lo) & (alpha*c-beta*lo <= 0), 0., q)
    return np.where(B >= cap, 0., q)


def _event_time(start, target, equilibrium, rate):
    ratio = np.divide(target - equilibrium, start - equilibrium,
                      out=np.ones_like(start), where=start != equilibrium)
    valid = (ratio > 0) & (ratio < 1) & (rate > 0)
    out = np.full_like(start, np.inf)
    out[valid] = -np.log(ratio[valid]) / rate[valid]
    return out


def _smooth_integrate(b, g, beta, lo, T, cap, delta, duration, p):
    """Adaptive SSP-RK3 with independent local error control for each value.

    The source coefficients and concentration remain frozen, exactly as in the
    sharp reference map. The h*L <= 1/2 safeguard limits the derivative of this
    smooth ODE; accuracy is separately controlled by step doubling (order 3).
    A physical B1/cap event completes a net source trajectory at that boundary.
    It never clips mobile inventory; the caller transfers delta B in a pair.
    """
    b = b.copy()
    rem = duration.copy()
    L = T - delta
    # |dq/dB| bound over the band and its adjacent analytic branches.
    lipschitz = (beta + g * (.75 / delta * (cap / L - 1) + cap / L**2)) / p.m0
    h = np.minimum(rem, np.divide(.5, lipschitz,
                                  out=np.full_like(rem, np.inf), where=lipschitz > 0))
    count = np.zeros(len(b), dtype=int)

    def rhs(v, gg, bb, tt, cc, dd):
        w = smooth_weight(v, tt, dd)
        quotient = np.divide(cc, v, out=np.zeros_like(v), where=v > 0)
        return (gg * (1 + w * (quotient - 1)) - bb * v) / p.m0

    def rk3(v, hh, gg, bb, tt, cc, dd):
        v1 = v + hh * rhs(v, gg, bb, tt, cc, dd)
        v2 = .75 * v + .25 * (v1 + hh * rhs(v1, gg, bb, tt, cc, dd))
        return (v + 2 * (v2 + hh * rhs(v2, gg, bb, tt, cc, dd))) / 3

    for _ in range(2000):
        active = rem > 0
        if not active.any():
            break
        ids = np.flatnonzero(active)
        hh = np.minimum(h[ids], rem[ids])
        args = (g[ids], beta[ids], T[ids], cap[ids], delta[ids])
        old = b[ids]
        full = rk3(old, hh, *args)
        half = rk3(old, hh / 2, *args)
        fine = rk3(half, hh / 2, *args)
        error = np.abs(fine - full) / 7
        tolerance = p.source_atol + p.source_rtol * np.maximum(np.abs(old), np.abs(fine))
        accept = (error <= tolerance) & np.isfinite(fine)
        if accept.any():
            jj = ids[accept]
            b[jj] = fine[accept]
            before = rem[jj].copy()
            rem[jj] = np.where(hh[accept] >= before, 0., before - hh[accept])
            count[jj] += 1
            DIAGNOSTICS['accepted_substeps'] += int(accept.sum())
            DIAGNOSTICS['maximum_accepted_estimated_error'] = max(
                DIAGNOSTICS['maximum_accepted_estimated_error'], float(error[accept].max()))
        DIAGNOSTICS['rejected_substeps'] += int((~accept).sum())
        factor = np.clip(.9 * np.power(np.divide(tolerance, error,
                              out=np.full_like(error, 1e8), where=error > 0), .25), .15, 3.)
        h[ids] = hh * factor
        h = np.minimum(h, np.divide(.5, lipschitz,
                                    out=np.full_like(h, np.inf), where=lipschitz > 0))
        if np.any((rem > 0) & (h <= np.finfo(float).tiny)):
            raise RuntimeError('Smooth source underflow in adaptive time step')
    else:
        raise RuntimeError('Smooth source adaptive integration did not finish')
    DIAGNOSTICS['maximum_substeps_per_value'] = max(
        DIAGNOSTICS['maximum_substeps_per_value'], int(count.max(initial=0)))
    DIAGNOSTICS['minimum_internal_deposit'] = min(
        DIAGNOSTICS['minimum_internal_deposit'], float(b.min(initial=np.inf)))
    lower = b < lo
    upper = b >= cap
    if np.any(lower & (g - beta * lo > 0)):
        raise RuntimeError('Adaptive map crossed an inward B1 boundary')
    if np.any(upper & (g - beta * cap < 0)):
        raise RuntimeError('Adaptive map crossed an inward capacity boundary')
    # Event completion of the scalar deposit ODE. The paired mobile update is
    # made once, using precisely this deposit increment, by the original solver.
    b[lower] = lo[lower]
    b[upper] = cap[upper]
    return b, int(lower.sum()), int(upper.sum())


def smooth_source_advance(B, c, alpha, beta, dt, p):
    DIAGNOSTICS['calls'] += 1
    g = np.broadcast_to(alpha * c, B.shape)
    beta = np.broadcast_to(beta, B.shape)
    lo = np.broadcast_to(np.asarray(p.B1)[:, None], B.shape)
    T = np.broadcast_to(np.asarray(p.B2)[:, None], B.shape)
    cap = np.broadcast_to(np.asarray(p.cap)[:, None], B.shape)
    chi = np.broadcast_to(np.asarray(p.charging)[:, None], B.shape)
    delta = T * p.smoothing_fraction
    L = T - delta
    H = T + delta
    rate = beta / p.m0
    eq = np.divide(g, beta, out=np.zeros_like(g), where=beta > 0)
    b = B.copy()
    rem = np.full_like(B, dt)
    holds = 0
    caps = 0
    for _ in range(6):
        active = rem > 0
        if not active.any():
            break
        rem[active & (b >= cap)] = 0.
        hold = active & (b == lo) & (g - beta * lo <= 0)
        holds += int(hold.sum())
        rem[hold] = 0.
        active = rem > 0
        charging = active & (b < lo)
        if charging.any():
            v = chi * g / p.m0
            event = np.divide(lo - b, v, out=np.full_like(b, np.inf), where=v > 0)
            step = np.minimum(rem, event)
            crossed = charging & (event <= rem)
            b[charging] += v[charging] * step[charging]
            b[crossed] = lo[crossed]
            rem[charging] -= step[charging]
        active = rem > 0
        qL = g - beta * L
        linear = active & ((b < L) | ((b == L) & (qL <= 0))) & (b >= lo)
        if linear.any():
            target = np.where(g - beta * b >= 0, L, lo)
            event = _event_time(b, target, eq, rate)
            zero = (beta == 0) & (g > 0)
            event[zero] = (target[zero] - b[zero]) * p.m0 / g[zero]
            step = np.minimum(rem, event)
            crossed = linear & (event <= rem)
            out = np.where(beta > 0, eq + (b - eq) * np.exp(-rate * step), b + g * step / p.m0)
            b[linear] = out[linear]
            b[crossed] = target[crossed]
            rem[linear] -= step[linear]
        active = rem > 0
        qH = g * cap / H - beta * H
        aging = active & ((b > H) | ((b == H) & (qH >= 0))) & (b < cap)
        if aging.any():
            target = np.where(g * cap - beta * b*b >= 0, cap, H)
            event = _event_time(b*b, target*target, eq*cap, 2*rate)
            zero = (beta == 0) & (g > 0)
            event[zero] = (target[zero]**2 - b[zero]**2) * p.m0 / (2*g[zero]*cap[zero])
            step = np.minimum(rem, event)
            crossed = aging & (event <= rem)
            out2 = np.where(beta > 0, eq*cap + (b*b-eq*cap)*np.exp(-2*rate*step),
                            b*b + 2*g*cap*step/p.m0)
            out = np.sqrt(np.maximum(0., out2))
            b[aging] = out[aging]
            b[crossed] = target[crossed]
            caps += int((crossed & (target == cap)).sum())
            rem[aging] -= step[aging]
        active = rem > 0
        inband = active & (b >= L) & (b <= H)
        if inband.any():
            DIAGNOSTICS['smooth_values'] += int(inband.sum())
            values, hc, cc = _smooth_integrate(
                b[inband], g[inband], beta[inband], lo[inband], T[inband], cap[inband],
                delta[inband], rem[inband], p)
            b[inband] = values
            rem[inband] = 0.
            holds += hc
            caps += cc
    if np.any(rem > 0):
        raise RuntimeError('Unresolved smooth source branch event')
    return b, holds, caps
