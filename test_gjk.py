import ctypes

import hypothesis.strategies as st
import numpy as np
from hypothesis import given


class Vec3(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_float),
        ("y", ctypes.c_float),
        ("z", ctypes.c_float),
    ]

    def __repr__(self):
        return f".{{.x={self.x}, .y={self.y}, .z={self.z}}}"


lib = ctypes.CDLL("zig-out/lib/libgjk.dylib")

vec3 = st.tuples(st.floats(), st.floats(), st.floats())


@st.composite
def vec3s(draw):
    x, y, z = draw(
        st.tuples(
            st.floats(min_value=-10000, max_value=10000),
            st.floats(min_value=-10000, max_value=10000),
            st.floats(min_value=-10000, max_value=10000),
        )
    )
    return Vec3(x, y, z)


@given(
    vec3s(),
    st.floats(min_value=0.01, max_value=1000),
    vec3s(),
    st.floats(min_value=0.01, max_value=1000),
)
def test_spheres(pos1, rad1, pos2, rad2):
    result = bool(
        lib.gjk_spheres(
            pos1,
            ctypes.c_float(rad1),
            pos2,
            ctypes.c_float(rad2),
        )
    )

    expected = (
        np.linalg.norm(
            np.array([pos1.x, pos1.y, pos1.z]) - np.array([pos2.x, pos2.y, pos2.z])
        )
        < rad1 + rad2
    )
    assert expected == result
