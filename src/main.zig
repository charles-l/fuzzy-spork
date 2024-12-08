const std = @import("std");
const BoundedArray = std.BoundedArray;
const print = std.debug.print;
const fmt = std.fmt;
const gjk = @import("gjk.zig");

const rl = @import("raylib");
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

var show_iter: usize = 11;
var debug_draw_count: usize = 0;

const Entity = struct {
    active: bool,

    position: rl.Vector3,
    scale: rl.Vector3,
    rotation: rl.Quaternion,
    collider: gjk.Collider,

    const Self = @This();
    fn fromPosition(position: rl.Vector3) Self {
        return std.mem.zeroInit(Self, .{
            .active = true,
            .position = position,
            .scale = .{ .x = 1, .y = 1, .z = 1 },
            .rotation = rl.Quaternion.identity(),
            .collider = .{ .sphere = .{ .radius = 0 } },
        });
    }
};

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

const max_entities = 16;
var entities = [_]Entity{std.mem.zeroInit(Entity, .{ .collider = .{ .sphere = .{ .radius = 0 } } })} ** max_entities;

const EntityIter = struct {
    i: usize,

    fn next(self: *@This()) ?*Entity {
        while (self.i < max_entities) {
            defer self.i += 1;
            if (entities[self.i].active) {
                return &entities[self.i];
            }
        }
        return null;
    }
};

fn entityIter() EntityIter {
    return .{ .i = 0 };
}

pub fn main() anyerror!void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer _ = gpa.deinit();

    var lua = try Lua.init(&allocator);
    defer lua.deinit();

    entities[0] = Entity.fromPosition(.{ .x = 0, .y = 0, .z = 0 });
    entities[0].collider = gjk.Collider{ .sphere = .{ .radius = 791 } };

    entities[1] = Entity.fromPosition(.{ .x = 0, .y = 3, .z = 791 });
    entities[1].collider = gjk.Collider{ .sphere = .{ .radius = 0.5 } };

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
        .position = entities[0].position.add(rl.Vector3.init(4, 2, 4)),
        .target = entities[0].position,
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

    const seed = 6;
    var rnd = std.rand.DefaultPrng.init(seed);
    const random_points = lbl: {
        var r: [20]vec3 = undefined;

        for (0..r.len) |i| {
            r[i] = gjk.vec3.init(
                (rnd.random().float(f32) - 0.5) * 5,
                (rnd.random().float(f32) - 0.5) * 5,
                (rnd.random().float(f32) - 0.5) * 5,
            );
        }

        break :lbl r;
    };

    //const cube_pos = rl.Vector3.init(0, 3, 0);

    //--------------------------------------------------------------------------------------

    rl.disableCursor();
    // Main game loop
    while (!rl.windowShouldClose()) { // Detect window close button or ESC key
        camera.update(rl.CameraMode.camera_free);

        if (rl.isKeyPressed(.key_z)) camera.target = entities[0].position;

        if (rl.isKeyDown(.key_l)) {
            entities[1].position.z += 0.1;
        }
        if (rl.isKeyDown(.key_h)) {
            entities[1].position.z -= 0.1;
        }
        if (rl.isKeyDown(.key_j)) {
            entities[1].position.y += 0.1;
        }
        if (rl.isKeyDown(.key_k)) {
            entities[1].position.y -= 0.1;
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

            for (random_points) |p| {
                rl.drawSphere(@bitCast(p), 0.1, rl.Color.dark_blue);
            }

            const maybe_hull = try gjk.quickhull(allocator, &random_points);
            if (maybe_hull) |hull| {
                allocator.free(hull.verts);
                allocator.free(hull.faces);
            }

            @import("debug.zig").debugDraw();

            if (false) {
                var iter1 = entityIter();
                while (iter1.next()) |entity1| {
                    var iter2 = EntityIter{ .i = iter1.i };
                    while (iter2.next()) |entity2| {
                        const r = gjk.gjk_inner(@bitCast(entity1.position), entity1.collider, @bitCast(entity2.position), entity2.collider, show_iter);
                        debug_draw_simplex(r.simplex);
                        debug_log("iters={}\n", .{show_iter});
                        if (r.collision) {
                            //const er = gjk.epa(r.simplex, @bitCast(entity1.position), entity1.collider, @bitCast(entity2.position), entity2.collider);
                            //print("{}\n", .{er});
                            debug_log("collisions\n", .{});
                            //debug_draw_simplex(r.simplex);
                        }
                    }
                }
            }

            if (false) {
                var it = entityIter();
                while (it.next()) |e| {
                    switch (e.collider) {
                        gjk.Collider.sphere => |s| {
                            rl.drawSphere(e.position, s.radius, rl.Color.purple);
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
