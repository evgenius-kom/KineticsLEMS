"""Multi-step kinetic topology enumeration.

The four ICTAC-recognised multi-step topologies plus a SINGLE marker for
the trivial one-step case.

When :func:`build_ode_system` is filled in, each branch should construct an
``f(state, T) → dstate/dt`` callable suitable for ``scipy.integrate.solve_ivp``:

* SINGLE         — state = [α], one channel.
* PARALLEL       — state = [α_1, …, α_n], independent; α_total = Σ w_i α_i.
* CONSECUTIVE    — state = species concentrations, A → B → C → … chain.
* COMPETITIVE    — shared reactant A → B and A → C, weights ≠ branching ratios.
* MIXED          — at least one parallel + one consecutive branch.
"""
from __future__ import annotations

from enum import StrEnum

import numpy as np

from .problem import FittingProblem


class Topology(StrEnum):
    SINGLE = "single"
    PARALLEL = "parallel"
    CONSECUTIVE = "consecutive"
    COMPETITIVE = "competitive"
    MIXED = "mixed"


def build_ode_system(problem: FittingProblem, params: np.ndarray):
    """Build the right-hand side function for the topology declared in ``problem``.

    Returns
    -------
    callable
        ``rhs(t, state, T_of_t)`` → ``dstate/dt`` for ``scipy.integrate.solve_ivp``.

    Raises
    ------
    NotImplementedError
        For every topology — this is infrastructure only. SINGLE will be the
        first wired up; PARALLEL needs the 2-parallel synthetic test bed (#16).
    """
    # TODO: implement SINGLE first.
    # Sketch for SINGLE:
    #
    #     def rhs(t, state, T_of_t):
    #         a = float(state[0])
    #         T = float(T_of_t(t))
    #         k = A * np.exp(-E / (R_GAS * T))
    #         return [k * f_alpha(a)]
    #
    # where (E, A) come from ``params`` and ``f_alpha`` from the reaction-model
    # registry in ``methods.master_plot.MASTER_MODELS``.
    raise NotImplementedError(
        f"Topology {problem.topology!r} ODE system is not yet implemented; "
        "see fitting/topology.py for the migration plan."
    )


__all__ = ["Topology", "build_ode_system"]
