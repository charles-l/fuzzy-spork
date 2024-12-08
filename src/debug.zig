const std = @import("std");
const rl = @import("raylib");

const max_faces = 1000;

var ctx = struct {
    step_target: usize = 1,
    steps: usize = 0,
    faces: std.BoundedArray([3]usize, max_faces) = .{ .len = 0 },
    verts: std.BoundedArray([3]f32, 1000) = .{ .len = 0 },
    highlight: std.StaticBitSet(max_faces) = std.StaticBitSet(max_faces).initEmpty(),
    // log file
    file: ?std.fs.File = null,
}{};

pub fn debugDraw() void {
    ctx.steps = 0;

    for (ctx.faces.slice(), 0..) |face, face_i| {
        const p1: rl.Vector3 = @bitCast(ctx.verts.slice()[face[0]]);
        const p2: rl.Vector3 = @bitCast(ctx.verts.slice()[face[1]]);
        const p3: rl.Vector3 = @bitCast(ctx.verts.slice()[face[2]]);

        const center = p1.add(p2).add(p3).mul(1.0 / 3.0);
        const norm = p3.sub(p1).cross(p2.sub(p1)).normalize();
        rl.drawLine3D(center, center.add(norm.mul(0.4)), rl.Color.green);

        rl.drawLine3D(@bitCast(p1), @bitCast(p2), rl.Color.blue);
        rl.drawLine3D(@bitCast(p2), @bitCast(p3), rl.Color.blue);
        rl.drawLine3D(@bitCast(p3), @bitCast(p1), rl.Color.blue);

        if (ctx.highlight.isSet(face_i)) {
            rl.drawTriangle3D(p3, p2, p1, rl.fade(rl.Color.yellow, 0.1));
        } else {
            rl.drawTriangle3D(p3, p2, p1, rl.fade(rl.Color.light_gray, 0.1));
        }
    }

    if (rl.isKeyPressed(rl.KeyboardKey.key_left_bracket) and ctx.step_target > 0) {
        if (rl.isKeyDown(rl.KeyboardKey.key_left_shift)) {
            ctx.step_target = 0;
        } else {
            ctx.step_target -= 1;
        }
    }

    if (rl.isKeyPressed(rl.KeyboardKey.key_right_bracket)) {
        ctx.step_target += 1;
    }
}

fn log(comptime T: type, writer: std.fs.File.Writer, prefix: []const u8, data: []const T) !void {
    _ = try writer.write(prefix);
    try writer.writeByte(':');
    try writer.writeInt(u32, @intCast(@sizeOf(T)), .little);
    try writer.writeInt(u32, @intCast(data.len), .little);
    for (data) |d| {
        const bytes: [@sizeOf(T)]u8 = @bitCast(d);
        _ = try writer.write(&bytes);
    }
}

pub fn debugMesh(verts: []const [3]f32, faces: []const [3]usize) void {
    if (ctx.file == null) {
        ctx.file = std.fs.createFileAbsolute("/tmp/log.txt", .{ .truncate = true, .read = true }) catch @panic("couldn't create log file");
    }
    const writer = ctx.file.?.writer();
    log([3]f32, writer, "verts", verts) catch @panic("log verts");
    log([3]usize, writer, "faces", faces) catch @panic("log faces");

    ctx.steps += 1;
    if (!isDebugIter()) {
        return;
    }

    ctx.faces.resize(0) catch unreachable;
    ctx.verts.resize(0) catch unreachable;
    ctx.highlight = std.StaticBitSet(max_faces).initEmpty();
    ctx.verts.appendSlice(verts) catch @panic("not enough space for verts");
    ctx.faces.appendSlice(faces) catch @panic("not enough space for faces");
}

pub fn isDebugIter() bool {
    return ctx.steps - 1 == ctx.step_target;
}

pub fn highlightFace(i: usize) void {
    if (isDebugIter()) {
        ctx.highlight.set(i);
    }
}
