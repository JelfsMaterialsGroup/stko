import numpy as np

from stko._internal.optimizers.gfnff import lbfgs_minimize


def test_lbfgs_minimize_quadratic() -> None:
    # Anisotropic quadratic with a known minimum: E = 1/2 sum k (x - x0)².
    rng = np.random.default_rng(4)
    minimum = rng.normal(size=(5, 3))
    k = rng.uniform(1.0, 50.0, size=(5, 3))

    def forces(positions: np.ndarray) -> np.ndarray:
        return -k * (positions - minimum)

    start = minimum + rng.normal(scale=2.0, size=(5, 3))
    positions, steps, largest_force = lbfgs_minimize(
        forces,
        start,
        fmax=1e-6,
        max_steps=500,
        max_step=0.2,
    )

    assert largest_force < 1e-6  # noqa: PLR2004
    assert steps < 500  # noqa: PLR2004
    assert np.allclose(positions, minimum, atol=1e-6)


def test_lbfgs_minimize_max_step() -> None:
    # From far away, no atom may move more than max_step per step.
    moves = []

    def forces(positions: np.ndarray) -> np.ndarray:
        moves.append(positions.copy())
        return -positions

    lbfgs_minimize(
        forces,
        np.full((3, 3), 10.0),
        fmax=1e-3,
        max_steps=200,
        max_step=0.2,
    )

    steps = np.diff(np.array(moves), axis=0)
    assert np.linalg.norm(steps, axis=2).max() <= 0.2 + 1e-12
