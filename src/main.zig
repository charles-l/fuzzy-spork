const rl = @import("raylib");
const BoundedArray = @import("std").BoundedArray;
const print = @import("std").debug.print;
const fmt = @import("std").fmt;

const vec3 = rl.Vector3;

const FrameLog = struct {
    buf: [1024]u8,
    offset: usize,
};

var frame_log = FrameLog{ .buf = undefined, .offset = 0 };

fn clearFrameLog() void {
    frame_log.offset = 0;
}

fn debugLog(comptime fmtstr: []const u8, args: anytype) void {
    const s = fmt.bufPrint(frame_log.buf[frame_log.offset..], fmtstr, args) catch "";
    frame_log.offset += s.len;
}

var show_iter: usize = 1;
var debug_draw_count: usize = 0;

pub fn main() anyerror!void {
    // Initialization
    //--------------------------------------------------------------------------------------
    const screenWidth = 800;
    const screenHeight = 450;

    rl.initWindow(screenWidth, screenHeight, "raylib-zig [core] example - basic window");
    defer rl.closeWindow(); // Close window and OpenGL context

    rl.setWindowFocused();

    rl.setTargetFPS(60); // Set our game to run at 60 frames-per-second

    var camera = rl.Camera3D{
        .position = rl.Vector3.init(4, 2, 4),
        .target = rl.Vector3.init(0, 0, 0),
        .up = rl.Vector3.init(0, 1, 0),
        .fovy = 60,
        .projection = rl.CameraProjection.camera_perspective,
    };

    const Sphere = struct {
        pos: rl.Vector3,
        radius: f32,
    };

    var spheres = [_]Sphere{
    //.{
    //    .pos = rl.Vector3.zero(),
    //    .radius = 1,
    //},
    .{
        .pos = rl.Vector3.init(0, 0, 3),
        .radius = 1,
    }};

    const cube_pos = rl.Vector3.init(0, 3, 0);

    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!rl.windowShouldClose()) { // Detect window close button or ESC key
        camera.update(rl.CameraMode.camera_free);

        if (rl.isKeyPressed(.key_z)) camera.target = rl.Vector3{ .x = 0, .y = 0, .z = 0 };

        if (rl.isKeyDown(.key_l)) {
            spheres[0].pos.z += 0.1;
        }
        if (rl.isKeyDown(.key_h)) {
            spheres[0].pos.z -= 0.1;
        }
        if (rl.isKeyDown(.key_j)) {
            spheres[0].pos.y += 0.1;
        }
        if (rl.isKeyDown(.key_k)) {
            spheres[0].pos.y -= 0.1;
        }

        if (rl.isKeyReleased(.key_left_bracket)) {
            if (show_iter > 1) {
                show_iter -= 1;
            }
        }
        if (rl.isKeyReleased(.key_right_bracket)) {
            show_iter += 1;
        }

        rl.beginDrawing();
        defer rl.endDrawing();
        defer clearFrameLog();
        defer {
            gjk_points.resize(0) catch unreachable;
            print("{}\n", .{debug_draw_count});
            debug_draw_count = 0;
        }

        rl.clearBackground(rl.Color.white);

        {
            camera.begin();
            defer camera.end();

            const r = gjk(
                Collider{
                    .pos = cube_pos,
                    .shape = .{ .mesh = .{ .verts = &.{
                        vec3.init(1, 1, 1),
                        vec3.init(1, 1, -1),
                        vec3.init(1, -1, 1),
                        vec3.init(1, -1, -1),
                        vec3.init(-1, 1, 1),
                        vec3.init(-1, 1, -1),
                        vec3.init(-1, -1, 1),
                        vec3.init(-1, -1, -1),
                    } } },
                },
                Collider{
                    .pos = spheres[0].pos,
                    .shape = .{ .sphere = .{ .radius = spheres[0].radius } },
                },
            );
            print("{}\n", .{r});

            rl.drawGrid(10, 1);

            rl.drawCube(cube_pos, 2, 2, 2, rl.Color.gray);
            for (spheres) |s| {
                rl.drawSphere(s.pos, s.radius, rl.Color.red.fade(if (r) 0.3 else 1));
            }
        }

        frame_log.buf[frame_log.offset] = 0;
        rl.drawText(frame_log.buf[0..frame_log.offset :0], 10, 10, 20, rl.Color.gray);
    }
}

var gjk_points = BoundedArray(vec3, 64).init(0) catch unreachable;
fn gjk(collider_a: Collider, collider_b: Collider) bool {
    // ref: https://winter.dev/articles/gjk-algorithm
    // ref: https://gist.github.com/vurtun/29727217c269a2fbf4c0ed9a1d11cb40
    var v = vec3.init(1, 0, 0);
    var a = collider_a.support(v).sub(collider_b.support(v.mul(-1)));
    var simplex = Simplex.new(&.{a});
    var d = a.mul(-1).normalize();

    gjk_points.append(a) catch {};

    rl.drawSphere(collider_a.support(v), 0.1, rl.Color.green);
    rl.drawSphere(collider_b.support(v.mul(-1)), 0.1, rl.Color.green);

    var i: usize = 0;

    simplex.debug_draw();

    for (0..10) |_| {
        i += 1;
        a = collider_a.support(d).sub(collider_b.support(d.mul(-1)));

        gjk_points.append(a) catch {};

        rl.drawSphere(collider_a.support(d.normalize()), 0.1, rl.Color.green);
        rl.drawSphere(collider_b.support(d.mul(-1).normalize()), 0.1, rl.Color.green);
        if (a.dot(d) < 0) {
            return false;
        }
        simplex.add(a);
        const contains_origin = next_simplex(&simplex, &d);
        d = d.normalize();
        if (contains_origin) {
            return true;
        }
    }

    debugLog("iters: {}\n", .{i});

    //print("shouldn't get here", .{});
    return false;
}

const ColliderType = enum {
    sphere,
    mesh,
};

const Collider = struct {
    pos: vec3,
    shape: union(ColliderType) {
        sphere: struct { radius: f32 },
        mesh: struct { verts: []const vec3 },
    },

    const Self = @This();
    fn support(self: Self, dir: vec3) vec3 {
        const point = switch (self.shape) {
            ColliderType.sphere => |s| dir.mul(s.radius),
            ColliderType.mesh => |m| lbl: {
                var max = m.verts[0];
                var max_dot = m.verts[0].dot(dir);
                for (m.verts[1..]) |v| {
                    const dot = v.dot(dir);
                    if (dot > max_dot) {
                        max = v;
                        max_dot = dot;
                    }
                }
                break :lbl max;
            },
        };
        return self.pos.add(point);
    }
};

const Simplex = struct {
    points: BoundedArray(vec3, 4),

    const Self = @This();

    fn new(vs: []const vec3) Self {
        return Self{ .points = BoundedArray(vec3, 4).fromSlice(vs) catch unreachable };
    }

    fn add(self: *Self, v: vec3) void {
        self.points.insert(0, v) catch @panic("too many points");
    }

    fn debug_draw(self: Self) void {
        debug_draw_count += 1;
        if (debug_draw_count == show_iter) {
            for (self.points.slice()) |p| {
                rl.drawSphere(p, 0.05, rl.Color.blue);
            }
            for (0..self.points.len) |i| {
                for (0..self.points.len) |j| {
                    if (i == j) {
                        continue;
                    }

                    rl.drawLine3D(self.points.get(i), self.points.get(j), rl.Color.blue);
                }
            }
        }
    }
};

fn next_simplex(simplex: *Simplex, dir: *vec3) bool {
    simplex.debug_draw();
    switch (simplex.points.len) {
        2 => {
            debugLog("line\n", .{});
            const a = simplex.points.get(0);
            const b = simplex.points.get(1);

            const ab = b.sub(a);
            const ao = a.mul(-1);
            if (ab.dot(ao) > 0) {
                dir.* = ab.cross(ao).cross(ab);
            } else {
                simplex.* = Simplex.new(&.{a});
                dir.* = ao;
            }

            return false;
        },
        3 => {
            debugLog("tri\n", .{});
            const a = simplex.points.get(0);
            const b = simplex.points.get(1);
            const c = simplex.points.get(2);

            const ab = b.sub(a);
            const ac = c.sub(a);
            const ao = a.mul(-1);

            const abc = ab.cross(ac);
            if (abc.cross(ac).dot(ao) > 0) {
                if (ac.dot(ao) > 0) {
                    simplex.* = Simplex.new(&.{ a, c });
                    dir.* = ac.cross(ao).cross(ac);
                } else {
                    simplex.* = Simplex.new(&.{ a, b });
                    return next_simplex(simplex, dir);
                }
            } else {
                if (ab.cross(abc).dot(ao) > 0) {
                    simplex.* = Simplex.new(&.{ a, b });
                    return next_simplex(simplex, dir);
                } else {
                    if (abc.dot(ao) > 0) {
                        dir.* = abc;
                    } else {
                        simplex.* = Simplex.new(&.{ a, c, b });
                        dir.* = abc.mul(-1);
                    }
                }
            }
            return false;
        },
        4 => {
            const a = simplex.points.get(0);
            const b = simplex.points.get(1);
            const c = simplex.points.get(2);
            const d = simplex.points.get(3);
            debugLog("simplex\n {}\n {}\n {}\n {}\n", .{ a, b, c, d });

            const ab = b.sub(a);
            const ac = c.sub(a);
            const ad = d.sub(a);
            const ao = a.mul(-1);

            const abc = ab.cross(ac);
            const acd = ac.cross(ad);
            const adb = ad.cross(ab);

            if (abc.dot(ao) > 0) {
                simplex.* = Simplex.new(&.{ a, b, c });
                return next_simplex(simplex, dir);
            }

            if (acd.dot(ao) > 0) {
                simplex.* = Simplex.new(&.{ a, c, d });
                return next_simplex(simplex, dir);
            }

            if (adb.dot(ao) > 0) {
                simplex.* = Simplex.new(&.{ a, d, b });
                return next_simplex(simplex, dir);
            }

            return true;
        },
        else => unreachable,
    }
}
