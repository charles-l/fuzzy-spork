import struct
from dataclasses import dataclass

import pyray as rl

# Initialization
SCREEN_WIDTH = 800
SCREEN_HEIGHT = 450

with open("log.txt", "rb") as f:
    log = f.read()


def eat_bytes(n):
    global log
    r = log[:n]
    log = log[n:]
    return r


parsed_log = []
# parse log
while True:
    try:
        tag, log = log.split(b":", 1)
        elem_size, length = struct.unpack("II", eat_bytes(8))
        data = eat_bytes(length * elem_size)
        if tag == b"verts":
            data = list(struct.iter_unpack("3f", data))
        elif tag == b"faces":
            data = list(struct.iter_unpack("3Q", data))
        else:
            assert False, f"unknown tag {tag}"
        assert len(data) == length, f"expected {length} got {len(data)}"
        print(tag, elem_size, length)
        parsed_log.append((tag.decode("ascii"), data))
    except ValueError:
        print("parsed")
        break

rl.init_window(SCREEN_WIDTH, SCREEN_HEIGHT, "raylib [core] example - 3d camera mode")

# Define the camera to look into our 3d world
camera = rl.Camera3D([0])
camera.position = rl.Vector3(0.0, 10.0, 10.0)  # Camera position
camera.target = rl.Vector3(0.0, 0.0, 0.0)  # Camera looking at point
camera.up = rl.Vector3(0.0, 1.0, 0.0)  # Camera up vector (rotation towards target)
camera.fovy = 45.0  # Camera field-of-view Y
camera.projection = rl.CameraProjection.CAMERA_PERSPECTIVE  # Camera mode type

cube_position = rl.Vector3(0.0, 0.0, 0.0)

rl.set_target_fps(60)

i = 0


@dataclass
class DebugObj:
    verts: list[tuple]
    faces: list[tuple]

    def draw(self):
        for vert in self.verts:
            rl.draw_sphere(vert, 0.05, rl.GRAY)
        for face in self.faces:
            p1, p2, p3 = (self.verts[i] for i in face)
            rl.draw_line_3d(p1, p2, rl.BLUE)
            rl.draw_line_3d(p2, p3, rl.BLUE)
            rl.draw_line_3d(p3, p1, rl.BLUE)
            # rl.draw_triangle_3d(p3, p2, p1, rl.fade(rl.GRAY, 0.1))


debug_obj = DebugObj([], [])
setattr(debug_obj, *parsed_log[0])

rl.disable_cursor()
# Main game loop
while not rl.window_should_close():
    # Draw
    rl.begin_drawing()

    rl.clear_background(rl.RAYWHITE)

    rl.begin_mode_3d(camera)
    rl.update_camera(camera, rl.CAMERA_FREE)

    rl.draw_grid(10, 1.0)
    debug_obj.draw()

    rl.end_mode_3d()

    rl.draw_fps(10, 10)

    rl.end_drawing()

    if rl.is_key_released(rl.KEY_LEFT_BRACKET):
        i = max(i - 1, 0)
        setattr(debug_obj, *parsed_log[i])
        print(parsed_log[i])
    elif rl.is_key_released(rl.KEY_RIGHT_BRACKET):
        i = min(i + 1, len(parsed_log) - 1)
        setattr(debug_obj, *parsed_log[i])
        print(parsed_log[i])


# De-Initialization
rl.close_window()
