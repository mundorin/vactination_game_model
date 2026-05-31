
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar


EPS = 1e-12


def _approx_equal(a: float, b: float, tol: float = 1e-8) -> bool:
    return abs(a - b) <= tol


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def herd_immunity_threshold(r0: float) -> float:
    if r0 <= 1.0:
        return 0.0
    return 1.0 - 1.0 / float(r0)


def final_size_sir_v(vaccination_rate: float, r0: float, vaccine_efficacy: float = 1.0) -> float:
    """Final epidemic size for a pre-epidemic SIR/V model in a well-mixed population."""
    x = _clamp01(vaccination_rate)
    susceptible_share = max(0.0, 1.0 - vaccine_efficacy * x)
    if susceptible_share <= EPS:
        return 0.0
    if r0 * susceptible_share <= 1.0 + 1e-12:
        return 0.0

    def equation(r: float) -> float:
        # use expm1 for numerical stability near r=0
        return r - susceptible_share * (-np.expm1(-r0 * r))

    left = 1e-8
    right = max(left * 10.0, susceptible_share)
    f_left = equation(left)
    f_right = equation(right)
    if f_left == 0.0:
        return float(left)
    if f_left * f_right > 0.0:
        # fallback to a simple grid search around the root in edge cases near the threshold
        grid = np.linspace(left, susceptible_share, 2000)
        vals = [equation(g) for g in grid]
        for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if va == 0.0:
                return float(a)
            if va * vb < 0.0:
                return float(brentq(equation, float(a), float(b), maxiter=200))
        return 0.0
    return float(brentq(equation, left, right, maxiter=200))


def infection_probability_unvaccinated(vaccination_rate: float, r0: float, vaccine_efficacy: float = 1.0) -> float:
    x = _clamp01(vaccination_rate)
    susceptible_share = max(0.0, 1.0 - vaccine_efficacy * x)
    if susceptible_share <= EPS:
        return 0.0
    return _clamp01(final_size_sir_v(x, r0=r0, vaccine_efficacy=vaccine_efficacy) / susceptible_share)


@dataclass
class SimulationConfig:
    max_time: float = 250.0
    n_eval: int = 1200
    rtol: float = 1e-8
    atol: float = 1e-10


class SIR_V:
    """SIR/V model for one isolated population with pre-epidemic vaccination."""

    def __init__(
        self,
        s: float,
        i: float,
        r: float,
        v: float,
        b: float,
        y: float,
        n: float = 1.0,
        treshhold: float = 1e-4,
        max_time: float = 250.0,
        n_eval: int = 1200,
    ) -> None:
        total = s + i + r + v
        if not _approx_equal(total, 1.0, tol=1e-5):
            raise ValueError(f"For SIR_V expected s+i+r+v=1, got {total:.6f}")
        if min(s, i, r, v) < -1e-12:
            raise ValueError("Compartments must be non-negative")
        if b < 0 or y <= 0:
            raise ValueError("Transmission and recovery rates must satisfy b>=0, y>0")

        self.s = float(s)
        self.i = float(i)
        self.r = float(r)
        self.v = float(v)
        self.b = float(b)
        self.y = float(y)
        self.n = float(n)
        self.treshhold = float(treshhold)
        self.config = SimulationConfig(max_time=max_time, n_eval=n_eval)
        self.time_step = max_time / max(1, n_eval - 1)
        self.timestamp = 0.0
        self.end = self.i <= self.treshhold
        self.model_info = self.make_simulation()

    def rhs(self, _t: float, z: np.ndarray) -> List[float]:
        s, i, r = z
        ds = -self.b * s * i
        di = self.b * s * i - self.y * i
        dr = self.y * i
        return [ds, di, dr]

    def make_simulation(self, max_time: Optional[float] = None, n_eval: Optional[int] = None) -> pd.DataFrame:
        max_time = self.config.max_time if max_time is None else float(max_time)
        n_eval = self.config.n_eval if n_eval is None else int(n_eval)
        t_eval = np.linspace(0.0, max_time, n_eval)
        sol = solve_ivp(
            self.rhs,
            (0.0, max_time),
            [self.s, self.i, self.r],
            t_eval=t_eval,
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        df = pd.DataFrame({
            "timestamp": sol.t,
            "s": sol.y[0],
            "i": sol.y[1],
            "r": sol.y[2],
            "v": np.full_like(sol.t, self.v),
        })
        self.timestamp = float(df["timestamp"].iloc[-1])
        self.end = float(df["i"].iloc[-1]) <= self.treshhold
        return df

    def iteretion(self) -> pd.Series:
        if self.timestamp >= self.config.max_time - EPS:
            self.end = True
            return self.model_info.iloc[-1]
        new_end = min(self.config.max_time, self.timestamp + self.time_step)
        t_eval = [self.timestamp, new_end]
        current = self.model_info.iloc[-1]
        sol = solve_ivp(
            self.rhs,
            (self.timestamp, new_end),
            [float(current["s"]), float(current["i"]), float(current["r"])],
            t_eval=t_eval,
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        row = {
            "timestamp": float(sol.t[-1]),
            "s": float(sol.y[0][-1]),
            "i": float(sol.y[1][-1]),
            "r": float(sol.y[2][-1]),
            "v": self.v,
        }
        self.model_info.loc[len(self.model_info)] = row
        self.timestamp = row["timestamp"]
        self.end = row["i"] <= self.treshhold or self.timestamp >= self.config.max_time - EPS
        return self.model_info.iloc[-1]

    def get_info(self) -> Dict[str, float]:
        last = self.model_info.iloc[-1]
        return {
            "s": float(last["s"]),
            "i": float(last["i"]),
            "r": float(last["r"]),
            "v": float(last["v"]),
            "b": self.b,
            "y": self.y,
            "n": self.n,
            "timestamp": float(last["timestamp"]),
            "end": bool(self.end),
        }

    def final_epidemic_size(self) -> float:
        return float(self.model_info["r"].iloc[-1])

    def infection_probability(self) -> float:
        hazard = self.b * self.model_info["i"].to_numpy()
        h = float(np.trapz(hazard, self.model_info["timestamp"].to_numpy()))
        return _clamp01(1.0 - math.exp(-h))

    def diplay_plots(self, title: str = "SIR/V dynamics") -> None:
        plt.figure(figsize=(7, 4.5))
        plt.plot(self.model_info["timestamp"], self.model_info["s"], label="S")
        plt.plot(self.model_info["timestamp"], self.model_info["i"], label="I")
        plt.plot(self.model_info["timestamp"], self.model_info["r"], label="R")
        plt.plot(self.model_info["timestamp"], self.model_info["v"], label="V")
        plt.xlabel("t")
        plt.ylabel("Fraction")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.show()


class SIRV:
    """SVIR model for one isolated population with dynamic vaccination."""

    def __init__(
        self,
        s: float,
        i: float,
        r: float,
        v: float,
        b: float,
        y: float,
        fi: float,
        n: float = 1.0,
        treshhold: float = 1e-4,
        max_time: float = 250.0,
        n_eval: int = 1200,
    ) -> None:
        total = s + i + r + v
        if not _approx_equal(total, 1.0, tol=1e-5):
            raise ValueError(f"For SIRV expected s+i+r+v=1, got {total:.6f}")
        if min(s, i, r, v) < -1e-12:
            raise ValueError("Compartments must be non-negative")
        if b < 0 or y <= 0 or fi < 0:
            raise ValueError("Rates must satisfy b>=0, y>0, fi>=0")

        self.s = float(s)
        self.i = float(i)
        self.r = float(r)
        self.v = float(v)
        self.b = float(b)
        self.y = float(y)
        self.fi = float(fi)
        self.n = float(n)
        self.treshhold = float(treshhold)
        self.config = SimulationConfig(max_time=max_time, n_eval=n_eval)
        self.time_step = max_time / max(1, n_eval - 1)
        self.timestamp = 0.0
        self.end = self.i <= self.treshhold
        self.model_info = self.make_simulation()

    def rhs(self, _t: float, z: np.ndarray) -> List[float]:
        s, v, i, r = z
        ds = -self.b * s * i - self.fi * s
        dv = self.fi * s
        di = self.b * s * i - self.y * i
        dr = self.y * i
        return [ds, dv, di, dr]

    def make_simulation(self, max_time: Optional[float] = None, n_eval: Optional[int] = None) -> pd.DataFrame:
        max_time = self.config.max_time if max_time is None else float(max_time)
        n_eval = self.config.n_eval if n_eval is None else int(n_eval)
        t_eval = np.linspace(0.0, max_time, n_eval)
        sol = solve_ivp(
            self.rhs,
            (0.0, max_time),
            [self.s, self.v, self.i, self.r],
            t_eval=t_eval,
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        df = pd.DataFrame({
            "timestamp": sol.t,
            "s": sol.y[0],
            "v": sol.y[1],
            "i": sol.y[2],
            "r": sol.y[3],
        })
        self.timestamp = float(df["timestamp"].iloc[-1])
        self.end = float(df["i"].iloc[-1]) <= self.treshhold
        return df

    def iteretion(self) -> pd.Series:
        if self.timestamp >= self.config.max_time - EPS:
            self.end = True
            return self.model_info.iloc[-1]
        new_end = min(self.config.max_time, self.timestamp + self.time_step)
        t_eval = [self.timestamp, new_end]
        current = self.model_info.iloc[-1]
        sol = solve_ivp(
            self.rhs,
            (self.timestamp, new_end),
            [float(current["s"]), float(current["v"]), float(current["i"]), float(current["r"])],
            t_eval=t_eval,
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        row = {
            "timestamp": float(sol.t[-1]),
            "s": float(sol.y[0][-1]),
            "v": float(sol.y[1][-1]),
            "i": float(sol.y[2][-1]),
            "r": float(sol.y[3][-1]),
        }
        self.model_info.loc[len(self.model_info)] = row
        self.timestamp = row["timestamp"]
        self.end = row["i"] <= self.treshhold or self.timestamp >= self.config.max_time - EPS
        return self.model_info.iloc[-1]

    def get_info(self) -> Dict[str, float]:
        last = self.model_info.iloc[-1]
        return {
            "s": float(last["s"]),
            "i": float(last["i"]),
            "r": float(last["r"]),
            "v": float(last["v"]),
            "b": self.b,
            "y": self.y,
            "fi": self.fi,
            "n": self.n,
            "timestamp": float(last["timestamp"]),
            "end": bool(self.end),
        }

    def diplay_plots(self, title: str = "SVIR dynamics") -> None:
        plt.figure(figsize=(7, 4.5))
        plt.plot(self.model_info["timestamp"], self.model_info["s"], label="S")
        plt.plot(self.model_info["timestamp"], self.model_info["i"], label="I")
        plt.plot(self.model_info["timestamp"], self.model_info["r"], label="R")
        plt.plot(self.model_info["timestamp"], self.model_info["v"], label="V")
        plt.xlabel("t")
        plt.ylabel("Fraction")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.show()


class SIRV_CHOICE(SIRV):
    """Compatibility wrapper around SVIR with the user's original class name."""
    pass


class SIR_V_CHOICE:
    """Single-city vaccination game with one-shot pre-epidemic vaccination."""

    def __init__(
        self,
        s: float,
        i: float,
        r: float,
        r_coef: float,
        c_i: float,
        c_v: float,
        treshhold: float = 1e-4,
        smoothness: float = 0.1,
        max_iter: int = 100000,
    ) -> None:
        total = s + i + r
        if not _approx_equal(total, 1.0, tol=1e-5):
            raise ValueError(f"For SIR_V_CHOICE expected s+i+r=1, got {total:.6f}")
        if min(s, i, r) < -1e-12:
            raise ValueError("Compartments must be non-negative")
        if r_coef <= 0 or c_i <= 0 or c_v < 0:
            raise ValueError("Need r_coef>0, c_i>0, c_v>=0")

        self.s = float(s)
        self.i = float(i)
        self.r = float(r)
        self.r_coef = float(r_coef)
        self.smoothness = float(smoothness)
        self.b = self.r_coef * self.smoothness
        self.y = self.smoothness
        self.c_i = float(c_i)
        self.c_v = float(c_v)
        self.treshhold = float(treshhold)
        self.max_iter = int(max_iter)
        self.time_step = 1.0
        self.timestamp = 0.0
        self.end = self.i <= self.treshhold

        self.v, self.p_infect, self.score, _, _ = self.get_nash_equilibrium()
        self.v_optimal, self.p_infect_optimal, self.score_optimal = self.get_pareto_equilibrium()
        self.model_info = self.make_simulation(v=self.v)

    @property
    def herd_threshold(self) -> float:
        return herd_immunity_threshold(self.r_coef)

    def get_infection_prob(self, v: float) -> float:
        return infection_probability_unvaccinated(vaccination_rate=v, r0=self.r_coef)

    def get_social_score(self, v: float, p_infect: float) -> float:
        v = _clamp01(v)
        p_infect = _clamp01(p_infect)
        return v * self.c_v + (1.0 - v) * p_infect * self.c_i

    def make_simulation(
        self,
        v: Optional[float] = None,
        max_time: Optional[float] = None,
        n_eval: int = 1200,
    ) -> pd.DataFrame:
        v = self.v if v is None else _clamp01(v)
        s0 = max(0.0, self.s - v)
        total = s0 + self.i + self.r + v
        if not _approx_equal(total, 1.0, tol=1e-5):
            if total > 0:
                s0 = s0 / total
                i0 = self.i / total
                r0 = self.r / total
                v0 = v / total
            else:
                raise ValueError("Degenerate initial state")
        else:
            i0 = self.i
            r0 = self.r
            v0 = v

        max_time = 250.0 if max_time is None else float(max_time)
        model = SIR_V(
            s=s0,
            i=i0,
            r=r0,
            v=v0,
            b=self.b,
            y=self.y,
            n=1.0,
            treshhold=self.treshhold,
            max_time=max_time,
            n_eval=n_eval,
        )
        return model.model_info.copy()

    def get_nash_equilibrium(self) -> Tuple[float, float, float, str, float]:
        target = self.c_v / self.c_i
        if target <= 0:
            nash_v = 1.0
            nash_p = 0.0
            score = self.get_social_score(nash_v, nash_p)
            return nash_v, nash_p, score, "full vaccination equilibrium", 0.0
        if target >= 1.0:
            nash_v = 0.0
            nash_p = self.get_infection_prob(0.0)
            score = self.get_social_score(nash_v, nash_p)
            return nash_v, nash_p, score, "no vaccination equilibrium", nash_p * self.c_i

        f0 = self.get_infection_prob(0.0) - target
        if f0 <= 0:
            nash_v = 0.0
            nash_p = self.get_infection_prob(nash_v)
            score = self.get_social_score(nash_v, nash_p)
            return nash_v, nash_p, score, "no vaccination equilibrium", nash_p * self.c_i

        upper = 0.999999
        f1 = self.get_infection_prob(upper) - target
        if f1 >= 0:
            nash_v = upper
            nash_p = self.get_infection_prob(nash_v)
            score = self.get_social_score(nash_v, nash_p)
            return nash_v, nash_p, score, "near-full vaccination equilibrium", nash_p * self.c_i

        root = brentq(lambda x: self.get_infection_prob(x) - target, 0.0, upper)
        nash_v = float(root)
        nash_p = self.get_infection_prob(nash_v)
        score = self.get_social_score(nash_v, nash_p)
        return nash_v, nash_p, score, "mixed equilibrium", nash_p * self.c_i

    def get_pareto_equilibrium(self) -> Tuple[float, float, float]:
        """Find the global minimum of the social-cost function C(x).

        With an ideal vaccine, C(x) is piecewise: it decreases on
        [0, 1 - 1/R0] (vaccinating cuts the epidemic) and increases on
        [1 - 1/R0, 1] (extra vaccinations beyond herd immunity are pure cost).
        This V-shape has a kink at x = 1 - 1/R0 and creates the *only* interior
        candidate minimum there.

        However, when c_v is large enough that herd-immunity vaccination costs
        more than the no-vaccination epidemic, the global minimum jumps to
        x = 0. So C(x) effectively has TWO local minima (at x = 0 and at the
        kink), and which one is global depends on c_v.

        scipy.optimize.minimize_scalar(method='bounded') gets stuck at the
        kink and misses the x = 0 alternative. Therefore we:
          1) evaluate C on a coarse grid,
          2) take the grid winner,
          3) refine locally around it with bounded scalar optimisation.

        This gives the correct global minimum across all values of c_v.
        """

        def cost(x: float) -> float:
            x = _clamp01(x)
            p = self.get_infection_prob(x)
            return self.get_social_score(x, p)

        # 1) Coarse grid scan. 200 points is enough to isolate which basin
        #    contains the global minimum (kink-basin or x = 0 basin).
        grid = np.linspace(0.0, 0.99, 200)
        costs = np.array([cost(float(x)) for x in grid])
        idx = int(np.argmin(costs))
        x_grid = float(grid[idx])

        # 2) Local refinement around the grid winner. Narrow window so
        #    minimize_scalar cannot escape into the other basin.
        lo = max(0.0, x_grid - 0.02)
        hi = min(0.999999, x_grid + 0.02)
        if hi > lo + 1e-9:
            result = minimize_scalar(
                lambda x: cost(float(x)),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-7},
            )
            x_refined = float(result.x)
            best_x = x_refined if cost(x_refined) < cost(x_grid) else x_grid
        else:
            best_x = x_grid

        best_p = self.get_infection_prob(best_x)
        best_score = self.get_social_score(best_x, best_p)
        return best_x, best_p, best_score

    def diplay_plots(self, title: str = "One-city equilibrium epidemic") -> None:
        plt.figure(figsize=(7, 4.5))
        plt.plot(self.model_info["timestamp"], self.model_info["s"], label="S")
        plt.plot(self.model_info["timestamp"], self.model_info["i"], label="I")
        plt.plot(self.model_info["timestamp"], self.model_info["r"], label="R")
        plt.plot(self.model_info["timestamp"], self.model_info["v"], label="V")
        plt.xlabel("t")
        plt.ylabel("Fraction")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def analyze_equilibria(self) -> Dict[str, Dict[str, float]]:
        pareto_v, pareto_p, pareto_cost = self.get_pareto_equilibrium()
        nash_v, nash_p, nash_cost, nash_condition, actual_condition = self.get_nash_equilibrium()
        efficiency_gap = nash_cost - pareto_cost
        is_efficient = abs(efficiency_gap) < 1e-8
        vaccination_gap = abs(nash_v - pareto_v)

        print("=" * 80)
        print("АНАЛИЗ РАВНОВЕСИЙ МОДЕЛИ ВАКЦИНАЦИИ")
        print("=" * 80)
        print("\nПАРАМЕТРЫ МОДЕЛИ:")
        print(f"R0: {self.r_coef:.3f}")
        print(f"gamma: {self.y:.3f}")
        print(f"beta: {self.b:.3f}")
        print(f"c_v: {self.c_v:.3f}")
        print(f"c_i: {self.c_i:.3f}")
        print(f"Порог коллективного иммунитета: {self.herd_threshold:.4f}")

        print("\nРАВНОВЕСИЕ НЭША:")
        print(f"v*: {nash_v:.4f}")
        print(f"P_infect(v*): {nash_p:.4f}")
        print(f"Социальные издержки при v*: {nash_cost:.4f}")
        print(f"Тип: {nash_condition}")
        print(f"Проверка условия c_i * P_infect = c_v: {actual_condition:.4f} vs {self.c_v:.4f}")

        print("\nСОЦИАЛЬНЫЙ ОПТИМУМ:")
        print(f"v_opt: {pareto_v:.4f}")
        print(f"P_infect(v_opt): {pareto_p:.4f}")
        print(f"Социальные издержки при v_opt: {pareto_cost:.4f}")

        print("\nСРАВНЕНИЕ:")
        print(f"Разница издержек: {efficiency_gap:.4f}")
        print(f"Разрыв по вакцинации: {vaccination_gap:.4f}")
        print(f"Эффективность: {'эффективно' if is_efficient else 'неэффективно'}")
        print("=" * 80)

        return {
            "parameters": {
                "R0": self.r_coef,
                "beta": self.b,
                "gamma": self.y,
                "c_v": self.c_v,
                "c_i": self.c_i,
                "herd_threshold": self.herd_threshold,
            },
            "nash_equilibrium": {
                "vaccination_rate": nash_v,
                "infection_prob": nash_p,
                "social_cost": nash_cost,
                "equilibrium_type": nash_condition,
                "condition_check": actual_condition,
            },
            "pareto_optimum": {
                "vaccination_rate": pareto_v,
                "infection_prob": pareto_p,
                "social_cost": pareto_cost,
            },
            "efficiency_analysis": {
                "cost_difference": efficiency_gap,
                "is_efficient": is_efficient,
                "vaccination_gap": vaccination_gap,
            },
        }


class TwoCityVaccinationGame:
    """Two-city SIR/V vaccination game from the diploma task."""

    def __init__(
        self,
        r_coef: float = 3.0,
        smoothness: float = 0.1,
        epsilon: float = 0.3,
        c_i: float = 1.0,
        c_v: float = 0.4,
        delta: float = 1e-4,
        max_time: float = 250.0,
        n_eval: int = 1200,
    ) -> None:
        if r_coef <= 0 or smoothness <= 0:
            raise ValueError("Need r_coef>0 and smoothness>0")
        if epsilon < 0:
            raise ValueError("epsilon must be non-negative")
        if c_i <= 0 or c_v < 0:
            raise ValueError("Need c_i>0 and c_v>=0")
        if not 0 <= delta < 1:
            raise ValueError("delta must lie in [0,1)")

        self.r_coef = float(r_coef)
        self.smoothness = float(smoothness)
        self.beta = self.r_coef * self.smoothness
        self.gamma = self.smoothness
        self.epsilon = float(epsilon)
        self.c_i = float(c_i)
        self.c_v = float(c_v)
        self.delta = float(delta)
        self.config = SimulationConfig(max_time=max_time, n_eval=n_eval)
        self._sim_cache: Dict[Tuple[float, float], pd.DataFrame] = {}

    def rhs(self, _t: float, z: np.ndarray) -> List[float]:
        s_a, i_a, r_a, s_b, i_b, r_b = z
        lambda_a = self.beta * (i_a + self.epsilon * i_b)
        lambda_b = self.beta * (i_b + self.epsilon * i_a)
        ds_a = -lambda_a * s_a
        di_a = lambda_a * s_a - self.gamma * i_a
        dr_a = self.gamma * i_a
        ds_b = -lambda_b * s_b
        di_b = lambda_b * s_b - self.gamma * i_b
        dr_b = self.gamma * i_b
        return [ds_a, di_a, dr_a, ds_b, di_b, dr_b]

    def initial_state(self, x: float, y: float) -> Tuple[float, float, float, float, float, float, float, float]:
        x = _clamp01(x)
        y = _clamp01(y)
        s_a0 = (1.0 - x) * (1.0 - self.delta)
        i_a0 = (1.0 - x) * self.delta
        r_a0 = 0.0
        v_a0 = x
        s_b0 = (1.0 - y) * (1.0 - self.delta)
        i_b0 = (1.0 - y) * self.delta
        r_b0 = 0.0
        v_b0 = y
        return s_a0, i_a0, r_a0, v_a0, s_b0, i_b0, r_b0, v_b0

    def simulate(self, x: float, y: float) -> pd.DataFrame:
        key = (round(float(x), 6), round(float(y), 6))
        if key in self._sim_cache:
            return self._sim_cache[key].copy()
        s_a0, i_a0, r_a0, v_a0, s_b0, i_b0, r_b0, v_b0 = self.initial_state(x, y)
        t_eval = np.linspace(0.0, self.config.max_time, self.config.n_eval)
        sol = solve_ivp(
            self.rhs,
            (0.0, self.config.max_time),
            [s_a0, i_a0, r_a0, s_b0, i_b0, r_b0],
            t_eval=t_eval,
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        df = pd.DataFrame({
            "timestamp": sol.t,
            "SA": sol.y[0],
            "IA": sol.y[1],
            "RA": sol.y[2],
            "VA": np.full_like(sol.t, v_a0),
            "SB": sol.y[3],
            "IB": sol.y[4],
            "RB": sol.y[5],
            "VB": np.full_like(sol.t, v_b0),
        })
        self._sim_cache[key] = df.copy()
        return df

    def infection_probability(self, x: float, y: float, city: str = "A") -> float:
        city = city.upper()
        df = self.simulate(x, y)
        if city == "A":
            hazard = self.beta * (df["IA"].to_numpy() + self.epsilon * df["IB"].to_numpy())
        elif city == "B":
            hazard = self.beta * (df["IB"].to_numpy() + self.epsilon * df["IA"].to_numpy())
        else:
            raise ValueError("city must be 'A' or 'B'")
        h = float(np.trapz(hazard, df["timestamp"].to_numpy()))
        return _clamp01(1.0 - math.exp(-h))

    def infection_probability_A(self, x: float, y: float) -> float:
        return self.infection_probability(x, y, city="A")

    def infection_probability_B(self, x: float, y: float) -> float:
        return self.infection_probability(x, y, city="B")

    def mean_social_cost_A(self, x: float, y: float) -> float:
        x = _clamp01(x)
        p = self.infection_probability_A(x, y)
        return x * self.c_v + (1.0 - x) * p * self.c_i

    def mean_social_cost_B(self, x: float, y: float) -> float:
        y = _clamp01(y)
        p = self.infection_probability_B(x, y)
        return y * self.c_v + (1.0 - y) * p * self.c_i

    def payoffs_A(self, x: float, y: float) -> Dict[str, float]:
        p = self.infection_probability_A(x, y)
        return {"Uvacc": -self.c_v, "Uno_vacc": -self.c_i * p}

    def nash_equilibrium(self, y: float) -> Tuple[float, str]:
        target = self.c_v / self.c_i
        if target <= 0:
            return 1.0, "full vaccination equilibrium"
        if target >= 1.0:
            return 0.0, "no vaccination equilibrium"

        f0 = self.infection_probability_A(0.0, y) - target
        if f0 <= 0:
            return 0.0, "boundary_0"

        upper = 0.999999
        f1 = self.infection_probability_A(upper, y) - target
        if f1 >= 0:
            return upper, "boundary_1"

        root = brentq(lambda x: self.infection_probability_A(x, y) - target, 0.0, upper)
        return float(root), "interior"

    def social_optimum(self, y: float) -> Tuple[float, float]:
        """Find the global minimum of mean_social_cost_A(x; y) over x ∈ [0, 1).

        Same robustness concern as SIR_V_CHOICE.get_pareto_equilibrium: with
        an ideal vaccine the cost function has a kink at the herd-immunity
        threshold and a competing minimum candidate at x = 0 when c_v is
        high. We do a grid scan first, then refine locally — this avoids the
        scipy `bounded` solver getting stuck on the wrong side of the kink.
        """

        def cost(x: float) -> float:
            x = max(0.0, min(0.999999, float(x)))
            return self.mean_social_cost_A(x, y)

        # 1) Coarse grid scan
        grid = np.linspace(0.0, 0.99, 200)
        costs = np.array([cost(float(x)) for x in grid])
        idx = int(np.argmin(costs))
        x_grid = float(grid[idx])

        # 2) Local refinement
        lo = max(0.0, x_grid - 0.02)
        hi = min(0.999999, x_grid + 0.02)
        if hi > lo + 1e-9:
            result = minimize_scalar(
                lambda x: cost(float(x)),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-7},
            )
            x_refined = float(result.x)
            x_opt = x_refined if cost(x_refined) < cost(x_grid) else x_grid
        else:
            x_opt = x_grid

        cost_opt = self.mean_social_cost_A(x_opt, y)
        return x_opt, cost_opt

    @staticmethod
    def normalized_sed(cost_nash: float, cost_opt: float) -> float:
        if cost_nash <= EPS:
            return 0.0
        return _clamp01((cost_nash - cost_opt) / cost_nash)

    def run_grid(self, y_grid: Iterable[float]) -> pd.DataFrame:
        rows: List[Dict[str, float]] = []
        for y in y_grid:
            y = _clamp01(float(y))
            x_nash, nash_type = self.nash_equilibrium(y)
            x_opt, cost_opt = self.social_optimum(y)
            p_nash = self.infection_probability_A(x_nash, y)
            p_opt = self.infection_probability_A(x_opt, y)
            cost_nash = self.mean_social_cost_A(x_nash, y)
            sed_raw = cost_nash - cost_opt
            sed_norm = self.normalized_sed(cost_nash, cost_opt)
            rows.append({
                "y": y,
                "x_nash": x_nash,
                "x_opt": x_opt,
                "p_nash": p_nash,
                "p_opt": p_opt,
                "cost_nash": cost_nash,
                "cost_opt": cost_opt,
                "sed_raw": sed_raw,
                "sed_norm": sed_norm,
                "nash_type": nash_type,
            })
        return pd.DataFrame(rows)

    def plot_example_dynamics(self, x: float, y: float, title: Optional[str] = None) -> None:
        df = self.simulate(x, y)
        title = title or f"Two-city dynamics: x={x:.2f}, y={y:.2f}"
        plt.figure(figsize=(8, 4.8))
        plt.plot(df["timestamp"], df["IA"], label="I_A")
        plt.plot(df["timestamp"], df["IB"], label="I_B")
        plt.plot(df["timestamp"], df["SA"], label="S_A", alpha=0.7)
        plt.plot(df["timestamp"], df["SB"], label="S_B", alpha=0.7)
        plt.xlabel("t")
        plt.ylabel("Fraction")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.show()


def simulate_sir_v_ode(
    s0: float,
    i0: float,
    r0_init: float,
    v: float,
    beta: float,
    gamma: float,
    t_max: float = 160.0,
    n_eval: int = 600,
) -> pd.DataFrame:
    model = SIR_V(s=s0, i=i0, r=r0_init, v=v, b=beta, y=gamma, max_time=t_max, n_eval=n_eval)
    return model.model_info.copy()


def simulate_svir_ode(
    s0: float,
    i0: float,
    r0_init: float,
    v0: float,
    phi: float,
    beta: float,
    gamma: float,
    t_max: float = 160.0,
    n_eval: int = 600,
) -> pd.DataFrame:
    model = SIRV(s=s0, i=i0, r=r0_init, v=v0, b=beta, y=gamma, fi=phi, max_time=t_max, n_eval=n_eval)
    return model.model_info.copy()


__all__ = [
    "SimulationConfig",
    "SIR_V",
    "SIRV",
    "SIRV_CHOICE",
    "SIR_V_CHOICE",
    "TwoCityVaccinationGame",
    "herd_immunity_threshold",
    "final_size_sir_v",
    "infection_probability_unvaccinated",
    "simulate_sir_v_ode",
    "simulate_svir_ode",
]
