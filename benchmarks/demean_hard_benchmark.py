"""
Hard FE demeaning benchmark suite.

Focus on problems that are difficult for iterative solvers:
- Conditioning problems (near-collinear, power-law, extreme imbalance)
- Structural problems (ladder/chain, hub/star, block diagonal)
- 3-FE scale problems (realistic firm-worker-year scenarios)

Usage:
    pixi run -e dev python benchmarks/demean_hard_benchmark.py
"""

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pyfixest.core._core_impl import _demean_rs

# =============================================================================
# Problem Generators
# =============================================================================


def near_collinear_2fe(
    n_obs: int,
    k: int = 500,
    correlation: float = 0.99,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create 2-FE problem where FE2 is nearly determined by FE1.

    This creates a near-singular Gram matrix, which is hard for GS.
    """
    rng = np.random.default_rng(seed)

    fe1 = rng.integers(0, k, size=n_obs)
    fe2 = fe1.copy()

    # Add noise: (1 - correlation) of observations get random FE2
    noise_pct = 1.0 - correlation
    noise_mask = rng.random(n_obs) < noise_pct
    fe2[noise_mask] = rng.integers(0, k, size=noise_mask.sum())

    flist = np.column_stack([fe1, fe2]).astype(np.uint64)

    # Create response
    true_alpha = rng.standard_normal(k) * 10
    true_beta = rng.standard_normal(k) * 5
    y = true_alpha[fe1] + true_beta[fe2] + rng.standard_normal(n_obs)

    desc = f"Near-collinear 2FE ({correlation * 100:.0f}% corr, k={k})"
    return flist, y.astype(np.float64), desc


def ladder_2fe(
    n_obs: int,
    k: int = 2000,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create bipartite ladder/chain structure.

    Each FE1 group connects to only 2 adjacent FE2 groups.
    This causes slow information propagation in coordinate descent.
    """
    rng = np.random.default_rng(seed)

    fe1 = np.arange(n_obs) % k
    # FE2 is fe1 or fe1+1 (mod k), creating a chain
    fe2 = (fe1 + rng.integers(0, 2, size=n_obs)) % k

    flist = np.column_stack([fe1, fe2]).astype(np.uint64)

    # Create response
    true_alpha = rng.standard_normal(k) * 10
    true_beta = rng.standard_normal(k) * 5
    y = true_alpha[fe1] + true_beta[fe2] + rng.standard_normal(n_obs)

    desc = f"Ladder 2FE (chain structure, k={k})"
    return flist, y.astype(np.float64), desc


def hub_2fe(
    n_obs: int,
    k: int = 1000,
    hub_fraction: float = 0.5,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create hub/star structure where one FE2 group dominates.

    Hub group (FE2=0) connects to all FE1 groups.
    Other FE2 groups connect to random FE1 groups.
    """
    rng = np.random.default_rng(seed)

    fe1 = rng.integers(0, k, size=n_obs)

    # Hub: FE2=0 gets hub_fraction of observations
    fe2 = np.zeros(n_obs, dtype=int)
    non_hub_mask = rng.random(n_obs) > hub_fraction
    fe2[non_hub_mask] = rng.integers(1, k, size=non_hub_mask.sum())

    flist = np.column_stack([fe1, fe2]).astype(np.uint64)

    # Create response
    true_alpha = rng.standard_normal(k) * 10
    true_beta = rng.standard_normal(k) * 5
    y = true_alpha[fe1] + true_beta[fe2] + rng.standard_normal(n_obs)

    desc = f"Hub 2FE ({hub_fraction * 100:.0f}% in hub, k={k})"
    return flist, y.astype(np.float64), desc


def power_law_2fe(
    n_obs: int,
    k: int = 1000,
    alpha: float = 1.5,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create 2-FE problem with power-law group size distribution.

    P(group i) ~ i^(-alpha)
    Creates high condition number due to unbalanced groups.
    """
    rng = np.random.default_rng(seed)

    # Power-law distribution for both FEs
    probs = np.power(np.arange(1, k + 1, dtype=float), -alpha)
    probs /= probs.sum()

    fe1 = rng.choice(k, size=n_obs, p=probs)
    fe2 = rng.choice(k, size=n_obs, p=probs)

    flist = np.column_stack([fe1, fe2]).astype(np.uint64)

    # Create response
    k1_actual = int(fe1.max()) + 1
    k2_actual = int(fe2.max()) + 1
    true_alpha = rng.standard_normal(k1_actual) * 10
    true_beta = rng.standard_normal(k2_actual) * 5
    y = true_alpha[fe1] + true_beta[fe2] + rng.standard_normal(n_obs)

    desc = f"Power-law 2FE (alpha={alpha}, k={k})"
    return flist, y.astype(np.float64), desc


def power_law_3fe(
    n_obs: int,
    k1: int = 500,
    k2: int = 500,
    k3: int = 100,
    alpha: float = 1.5,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Create 3-FE problem with power-law group sizes (firm-worker-year style)."""
    rng = np.random.default_rng(seed)

    def power_law_sample(k: int, n: int) -> np.ndarray:
        probs = np.power(np.arange(1, k + 1, dtype=float), -alpha)
        probs /= probs.sum()
        return rng.choice(k, size=n, p=probs)

    fe1 = power_law_sample(k1, n_obs)
    fe2 = power_law_sample(k2, n_obs)
    fe3 = power_law_sample(k3, n_obs)

    flist = np.column_stack([fe1, fe2, fe3]).astype(np.uint64)

    # Create response
    k1_actual = int(fe1.max()) + 1
    k2_actual = int(fe2.max()) + 1
    k3_actual = int(fe3.max()) + 1
    true_a = rng.standard_normal(k1_actual) * 10
    true_b = rng.standard_normal(k2_actual) * 5
    true_c = rng.standard_normal(k3_actual) * 2
    y = true_a[fe1] + true_b[fe2] + true_c[fe3] + rng.standard_normal(n_obs)

    desc = f"Power-law 3FE (alpha={alpha}, k={k1},{k2},{k3})"
    return flist, y.astype(np.float64), desc


def extreme_imbalance_3fe(
    n_obs: int,
    dominant_pct: float = 0.9,
    k1: int = 100,
    k2: int = 100,
    k3: int = 50,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create 3-FE problem where one group dominates in FE1.

    Group 0 gets dominant_pct of observations.
    """
    rng = np.random.default_rng(seed)

    # FE1: Group 0 gets dominant_pct, rest spread over k1-1 groups
    fe1 = np.zeros(n_obs, dtype=int)
    non_dominant = int(n_obs * (1 - dominant_pct))
    fe1[:non_dominant] = rng.integers(1, k1, size=non_dominant)

    fe2 = rng.integers(0, k2, size=n_obs)
    fe3 = rng.integers(0, k3, size=n_obs)

    flist = np.column_stack([fe1, fe2, fe3]).astype(np.uint64)

    # Create response
    true_a = rng.standard_normal(k1) * 10
    true_b = rng.standard_normal(k2) * 5
    true_c = rng.standard_normal(k3) * 2
    y = true_a[fe1] + true_b[fe2] + true_c[fe3] + rng.standard_normal(n_obs)

    desc = f"Extreme imbalance 3FE ({dominant_pct * 100:.0f}% dominant)"
    return flist, y.astype(np.float64), desc


def block_diagonal_3fe(
    n_obs: int,
    n_blocks: int = 10,
    k_per_block: int = 10,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Create block-diagonal structure (disconnected subproblems).

    Each block has its own set of FE groups.
    """
    rng = np.random.default_rng(seed)
    block_size = n_obs // n_blocks

    fe1 = np.zeros(n_obs, dtype=int)
    fe2 = np.zeros(n_obs, dtype=int)
    fe3 = np.zeros(n_obs, dtype=int)

    for b in range(n_blocks):
        start = b * block_size
        end = start + block_size if b < n_blocks - 1 else n_obs
        block_len = end - start

        base = b * k_per_block
        fe1[start:end] = base + rng.integers(0, k_per_block, size=block_len)
        fe2[start:end] = base + rng.integers(0, k_per_block, size=block_len)
        fe3[start:end] = base + rng.integers(0, k_per_block, size=block_len)

    flist = np.column_stack([fe1, fe2, fe3]).astype(np.uint64)

    # Create response
    n_coef = n_blocks * k_per_block
    true_a = rng.standard_normal(n_coef) * 10
    true_b = rng.standard_normal(n_coef) * 5
    true_c = rng.standard_normal(n_coef) * 2
    y = true_a[fe1] + true_b[fe2] + true_c[fe3] + rng.standard_normal(n_obs)

    desc = f"Block-diagonal 3FE ({n_blocks} blocks)"
    return flist, y.astype(np.float64), desc


def many_small_groups_3fe(
    n_obs: int,
    avg_obs_per_group: tuple[int, int, int] = (10, 20, 50),
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Create 3-FE problem with many small groups (high coefficient count)."""
    rng = np.random.default_rng(seed)

    k1 = n_obs // avg_obs_per_group[0]
    k2 = n_obs // avg_obs_per_group[1]
    k3 = n_obs // avg_obs_per_group[2]

    fe1 = rng.integers(0, k1, size=n_obs)
    fe2 = rng.integers(0, k2, size=n_obs)
    fe3 = rng.integers(0, k3, size=n_obs)

    flist = np.column_stack([fe1, fe2, fe3]).astype(np.uint64)

    # Create response
    true_a = rng.standard_normal(k1) * 10
    true_b = rng.standard_normal(k2) * 5
    true_c = rng.standard_normal(k3) * 2
    y = true_a[fe1] + true_b[fe2] + true_c[fe3] + rng.standard_normal(n_obs)

    n_coef = k1 + k2 + k3
    desc = f"Many small groups 3FE ({n_coef:,} coefficients)"
    return flist, y.astype(np.float64), desc


def dense_cross_3fe(
    n_obs: int,
    k1: int = 100,
    k2: int = 100,
    k3: int = 50,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Create fully crossed design where all (FE1, FE2) pairs are observed."""
    rng = np.random.default_rng(seed)

    n_pairs = k1 * k2
    obs_per_pair = max(1, n_obs // n_pairs)

    fe1 = np.tile(np.repeat(np.arange(k1), k2), obs_per_pair)[:n_obs]
    fe2 = np.tile(np.tile(np.arange(k2), k1), obs_per_pair)[:n_obs]

    # Pad if needed
    if len(fe1) < n_obs:
        extra = n_obs - len(fe1)
        fe1 = np.concatenate([fe1, rng.integers(0, k1, size=extra)])
        fe2 = np.concatenate([fe2, rng.integers(0, k2, size=extra)])

    fe3 = rng.integers(0, k3, size=n_obs)

    flist = np.column_stack([fe1, fe2, fe3]).astype(np.uint64)

    # Create response
    true_a = rng.standard_normal(k1) * 10
    true_b = rng.standard_normal(k2) * 5
    true_c = rng.standard_normal(k3) * 2
    y = true_a[fe1] + true_b[fe2] + true_c[fe3] + rng.standard_normal(n_obs)

    desc = f"Dense cross 3FE ({k1}x{k2} pairs)"
    return flist, y.astype(np.float64), desc


# =============================================================================
# Benchmark Runner
# =============================================================================

# Solvers to compare
SOLVERS = [
    ("gauss_seidel", "none", "GS"),
    ("lsmr", "none", "LSMR"),
    ("lsmr", "diagonal", "LSMR+diag"),
    ("lsmr", "streaming", "LSMR+stream"),
    ("lsmr", "sparse_gram", "LSMR+sparse"),
]


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    problem: str
    solver: str
    preconditioner: str
    label: str
    time_ms: float
    iterations: int
    converged: bool
    n_obs: int
    n_coef: int


def run_problem(
    flist: np.ndarray,
    y: np.ndarray,
    solver: str,
    preconditioner: str,
    tol: float = 1e-8,
    maxiter: int = 10000,
    n_runs: int = 3,
) -> tuple[float, int, bool]:
    """Run solver and return (time_ms, iterations, converged)."""
    x = y.reshape(-1, 1).astype(np.float64)

    times = []
    iterations = 0
    converged = False

    for _ in range(n_runs):
        start = time.perf_counter()
        result = _demean_rs(
            x=x,
            flist=flist,
            weights=None,
            tol=tol,
            maxiter=maxiter,
            reorder_fe=True,
            solver=solver,
            lsmr_preconditioner=preconditioner,
        )
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        iterations = result["iterations"]
        converged = result["success"]

    return min(times) * 1000, iterations, converged


def run_benchmark_suite(
    problems: list[tuple[np.ndarray, np.ndarray, str]],
    tol: float = 1e-8,
    maxiter: int = 10000,
) -> pd.DataFrame:
    """Run full benchmark suite."""
    results = []

    for flist, y, desc in problems:
        n_obs = len(y)
        n_coef = sum(flist[:, i].max() + 1 for i in range(flist.shape[1]))

        print(f"\n{desc}")
        print(f"  n_obs={n_obs:,}, n_coef={n_coef:,}")
        print("-" * 70)

        case_results = []
        for solver, precond, label in SOLVERS:
            try:
                time_ms, iters, converged = run_problem(
                    flist, y, solver, precond, tol, maxiter
                )
                result = BenchmarkResult(
                    problem=desc,
                    solver=solver,
                    preconditioner=precond,
                    label=label,
                    time_ms=time_ms,
                    iterations=iters,
                    converged=converged,
                    n_obs=n_obs,
                    n_coef=n_coef,
                )
                case_results.append(result)

                status = "ok" if converged else "FAIL"
                print(f"  {label:15s}: {time_ms:8.1f} ms, {iters:5d} iters [{status}]")

            except Exception as e:
                print(f"  {label:15s}: ERROR - {e}")

        results.extend(case_results)

        # Find winner
        converged_results = [r for r in case_results if r.converged]
        if converged_results:
            best = min(converged_results, key=lambda r: r.time_ms)
            gs = next((r for r in case_results if r.label == "GS"), None)
            if gs and gs.converged and best.label != "GS":
                speedup = gs.time_ms / best.time_ms
                print(f"  -> Winner: {best.label} ({speedup:.2f}x vs GS)")
            else:
                print(f"  -> Winner: {best.label}")

    return pd.DataFrame([r.__dict__ for r in results])


# =============================================================================
# Problem Sets
# =============================================================================


def get_conditioning_problems(n_obs: int = 200_000) -> list:
    """Problems that stress conditioning (GS should struggle)."""
    return [
        near_collinear_2fe(n_obs, k=500, correlation=0.99),
        near_collinear_2fe(n_obs, k=1000, correlation=0.995),
        power_law_2fe(n_obs, k=1000, alpha=1.5),
        power_law_2fe(n_obs, k=2000, alpha=2.0),
        extreme_imbalance_3fe(n_obs, dominant_pct=0.9),
        extreme_imbalance_3fe(n_obs, dominant_pct=0.95),
    ]


def get_structural_problems(n_obs: int = 200_000) -> list:
    """Problems with difficult graph structures."""
    return [
        ladder_2fe(n_obs, k=2000),
        ladder_2fe(n_obs, k=5000),
        hub_2fe(n_obs, k=1000, hub_fraction=0.5),
        hub_2fe(n_obs, k=1000, hub_fraction=0.8),
        block_diagonal_3fe(n_obs, n_blocks=10, k_per_block=10),
        block_diagonal_3fe(n_obs, n_blocks=20, k_per_block=20),
    ]


def get_scale_problems(n_obs: int = 500_000) -> list:
    """Large 3-FE problems (realistic firm-worker-year)."""
    return [
        power_law_3fe(n_obs, k1=500, k2=500, k3=100, alpha=1.5),
        power_law_3fe(n_obs, k1=1000, k2=1000, k3=200, alpha=1.5),
        many_small_groups_3fe(n_obs, avg_obs_per_group=(10, 20, 50)),
        dense_cross_3fe(n_obs, k1=100, k2=100, k3=50),
    ]


def get_all_problems() -> list:
    """Get all benchmark problems."""
    return (
        get_conditioning_problems(200_000)
        + get_structural_problems(200_000)
        + get_scale_problems(500_000)
    )


# =============================================================================
# Main
# =============================================================================


def print_summary(df: pd.DataFrame) -> None:
    """Print summary statistics."""
    print("\n" + "=" * 90)
    print("SUMMARY - Time (ms)")
    print("=" * 90)

    if df.empty:
        print("No results.")
        return

    pivot = df.pivot_table(
        values="time_ms", index="problem", columns="label", aggfunc="first"
    )
    col_order = ["GS", "LSMR", "LSMR+diag", "LSMR+stream", "LSMR+sparse"]
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    print(pivot.round(1).to_string())

    # Speedup vs GS
    if "GS" in pivot.columns:
        print("\n" + "=" * 90)
        print("SPEEDUP vs Gauss-Seidel (>1 = LSMR faster)")
        print("=" * 90)

        for col in pivot.columns:
            if col != "GS":
                speedup = pivot["GS"] / pivot[col]
                wins = (speedup > 1).sum()
                avg = speedup.mean()
                best = speedup.max()
                print(
                    f"{col}: avg {avg:.2f}x, best {best:.2f}x ({wins}/{len(speedup)} wins)"
                )

    # Iteration counts
    print("\n" + "=" * 90)
    print("ITERATIONS")
    print("=" * 90)

    pivot_iter = df.pivot_table(
        values="iterations", index="problem", columns="label", aggfunc="first"
    )
    pivot_iter = pivot_iter[[c for c in col_order if c in pivot_iter.columns]]
    print(pivot_iter.to_string())


def main():
    """Run the hard FE demeaning benchmark suite."""
    print("=" * 90)
    print("HARD FE DEMEANING BENCHMARK SUITE")
    print("=" * 90)
    print("\nFocusing on problems where GS may struggle:")
    print("  - Conditioning: near-collinear, power-law, extreme imbalance")
    print("  - Structure: ladder/chain, hub/star, block-diagonal")
    print("  - Scale: 500K obs, many coefficients")

    problems = get_all_problems()
    df = run_benchmark_suite(problems)

    print_summary(df)

    # Save results
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "demean_hard_benchmark.csv"
    df.to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")

    return df


if __name__ == "__main__":
    df = main()
