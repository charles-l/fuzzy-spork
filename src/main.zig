const std = @import("std");
const BoundedArray = std.BoundedArray;
const print = std.debug.print;
const fmt = std.fmt;
const gjk = @import("gjk.zig");

const rl = @import("raylib");
const ecs = @import("zigecs");
const ziglua = @import("ziglua");

const Lua = ziglua.Lua;

const vec3 = gjk.vec3;

const FrameLog = struct {
    buf: [1024]u8,
    offset: usize,
};

var frame_log = FrameLog{ .buf = undefined, .offset = 0 };

fn clearFrameLog() void {
    frame_log.offset = 0;
}

fn debug_log(comptime fmtstr: []const u8, args: anytype) void {
    const s = fmt.bufPrint(frame_log.buf[frame_log.offset..], fmtstr, args) catch "";
    frame_log.offset += s.len;
}

var show_iter: usize = 10;
var debug_draw_count: usize = 0;

const Position = rl.Vector3;

fn debug_draw_simplex(self: gjk.Simplex) void {
    for (self.points.slice()) |p| {
        rl.drawSphere(@bitCast(p), 0.05, rl.Color.blue);
    }
    for (0..self.points.len) |i| {
        for (0..self.points.len) |j| {
            if (i == j) {
                continue;
            }

            rl.drawLine3D(@bitCast(self.points.get(i)), @bitCast(self.points.get(j)), rl.Color.blue);
        }
    }
}

fn collision_system(positions: []Position, colliders: []gjk.Collider) void {
    for (positions, colliders, 0..) |pos1, col1, i| {
        for (positions[i + 1 ..], colliders[i + 1 ..]) |pos2, col2| {
            const r = gjk(pos1, col1, pos2, col2);
            if (r) {
                print("collisions\n", .{});
            }
        }
    }
}

pub fn main() anyerror!void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer _ = gpa.deinit();

    var lua = try Lua.init(&allocator);
    defer lua.deinit();

    var reg = ecs.Registry.init(allocator);
    defer reg.deinit();

    // pos1 = .{.x=32769.0, .y=0.0, .z=0.0}, rad1 = 0.001953125, pos2 = .{.x=32769.0, .y=0.0, .z=0.0}, rad2 = 0.001953125

    const spherePos = Position{ .x = 32769.0, .y = 0.0, .z = 0.0 };
    const sphere = reg.create();
    {
        reg.add(sphere, spherePos);
        reg.add(sphere, gjk.Collider{
            .sphere = .{ .radius = 0.1 },
        });
    }

    const sphere2 = reg.create();
    {
        reg.add(sphere2, Position{ .x = 32769.0, .y = 0.0, .z = 0.0 });
        reg.add(sphere2, gjk.Collider{
            .sphere = .{ .radius = 0.1 },
        });
    }

    //const cube = reg.create();
    //{
    //    reg.add(cube, Position.init(0, 1, 0));
    //    reg.add(cube, gjk.Collider{
    //        .mesh = .{ .verts = &.{
    //            vec3.init(1, 1, 1),
    //            vec3.init(1, 1, -1),
    //            vec3.init(1, -1, 1),
    //            vec3.init(1, -1, -1),
    //            vec3.init(-1, 1, 1),
    //            vec3.init(-1, 1, -1),
    //            vec3.init(-1, -1, 1),
    //            vec3.init(-1, -1, -1),
    //        } },
    //    });
    //}

    // Initialization
    //--------------------------------------------------------------------------------------
    const screenWidth = 800;
    const screenHeight = 450;

    rl.initWindow(screenWidth, screenHeight, "raylib-zig [core] example - basic window");
    defer rl.closeWindow(); // Close window and OpenGL context

    rl.setWindowFocused();

    rl.setTargetFPS(60); // Set our game to run at 60 frames-per-second

    var camera = rl.Camera3D{
        .position = spherePos.add(rl.Vector3.init(4, 2, 4)),
        .target = spherePos,
        .up = rl.Vector3.init(0, 1, 0),
        .fovy = 60,
        .projection = rl.CameraProjection.camera_perspective,
    };

    const origin_cam = rl.Camera3D{
        .position = rl.Vector3.init(4, 2, 4),
        .target = rl.Vector3.init(0, 0, 0),
        .up = rl.Vector3.init(0, 1, 0),
        .fovy = 60,
        .projection = rl.CameraProjection.camera_perspective,
    };

    //const cube_pos = rl.Vector3.init(0, 3, 0);

    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!rl.windowShouldClose()) { // Detect window close button or ESC key
        var view = reg.view(.{ Position, gjk.Collider }, .{});
        camera.update(rl.CameraMode.camera_free);

        if (rl.isKeyPressed(.key_z)) camera.target = spherePos;

        if (rl.isKeyDown(.key_l)) {
            view.get(Position, sphere).z += 0.1;
        }
        if (rl.isKeyDown(.key_h)) {
            view.get(Position, sphere).z -= 0.1;
        }
        if (rl.isKeyDown(.key_j)) {
            view.get(Position, sphere).y += 0.1;
        }
        if (rl.isKeyDown(.key_k)) {
            view.get(Position, sphere).y -= 0.1;
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

        rl.clearBackground(rl.Color.white);

        {
            const cam = if (rl.isKeyDown(.key_z)) origin_cam else camera;
            cam.begin();
            defer cam.end();

            rl.drawGrid(10, 1);

            {
                const positions = reg.raw(Position);
                const colliders = reg.raw(gjk.Collider);
                for (positions, colliders, 0..) |pos1, col1, i| {
                    for (positions[i + 1 ..], colliders[i + 1 ..]) |pos2, col2| {
                        const r = gjk.gjk_inner(@bitCast(pos1), col1, @bitCast(pos2), col2, show_iter);
                        if (r.collision) {
                            debug_log("collisions\n", .{});
                            debug_draw_simplex(r.simplex);
                        }
                    }
                }
            }

            //rl.drawCube(cube_pos, 2, 2, 2, rl.Color.gray);
            //for (spheres) |s| {
            //    //rl.drawSphere(s.pos, s.radius, rl.Color.red.fade(if (r) 0.3 else 1));
            //}
            //
            {
                var it = view.entityIterator();
                while (it.next()) |e| {
                    const position = view.getConst(Position, e);
                    switch (view.getConst(gjk.Collider, e)) {
                        gjk.Collider.sphere => |s| {
                            rl.drawSphere(position, s.radius, rl.Color.purple);
                        },
                        else => {},
                    }
                }
            }
        }

        frame_log.buf[frame_log.offset] = 0;
        rl.drawText(frame_log.buf[0..frame_log.offset :0], 10, 10, 20, rl.Color.gray);
    }
}
