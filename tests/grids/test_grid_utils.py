import numpy as np
import pytest
from layowt.grids.utils import _rotation_matrix


@pytest.mark.parametrize(
    ("angle"),
    np.arange(-360, 361, 1)
)
def test_rotation_matrix(angle):
    """
    Tests if the rotation matrix is correct (duplicating code, need to rethink).
    """
    # TODO: Need to rethink this test.
    result = [
            [np.cos(angle * np.pi / 180), -np.sin(angle * np.pi / 180)],
            [np.sin(angle * np.pi / 180), np.cos(angle * np.pi / 180)],
            ]

    np.testing.assert_almost_equal(_rotation_matrix(angle), np.array(result))

@pytest.mark.parametrize(
    ("angle"),
    np.arange(-360, 361, 1)
)
def test_rotation_matrix_orthogonal(angle):
    """
    Tests that the rotation matrix is orthogonal.
    """
    matrix = _rotation_matrix(angle)
    
    assert np.linalg.det(matrix) == pytest.approx(1, 0.0000001)
