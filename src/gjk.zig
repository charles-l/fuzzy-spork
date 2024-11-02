const std = @import("std");
pub const vec3 = extern struct {
    x: f32,
    y: f32,
    z: f32,

    const Self = @This();
    pub fn init(x: f32, y: f32, z: f32) Self {
        return .{ .x = x, .y = y, .z = z };
    }

    pub fn add(self: Self, v: vec3) Self {
        return .{ .x = self.x + v.x, .y = self.y + v.y, .z = self.z + v.z };
    }

    pub fn sub(self: Self, v: vec3) Self {
        return .{ .x = self.x - v.x, .y = self.y - v.y, .z = self.z - v.z };
    }

    pub fn mul(self: Self, f: f32) Self {
        return .{ .x = self.x * f, .y = self.y * f, .z = self.z * f };
    }

    pub fn dot(self: Self, v: vec3) f32 {
        return self.x * v.x + self.y * v.y + self.z * v.z;
    }

    pub fn normalize(self: Self) Self {
        const len = @sqrt(self.dot(self));
        return .{ .x = self.x / len, .y = self.y / len, .z = self.z / len };
    }

    pub fn cross(self: Self, v: vec3) Self {
        return .{
            .x = self.y * v.z - v.y * self.z,
            .y = self.z * v.x - v.z * self.x,
            .z = self.x * v.y - v.x * self.y,
        };
    }
};

const BoundedArray = @import("std").BoundedArray;

export fn gjk_spheres(a_pos: vec3, a_rad: f32, b_pos: vec3, b_rad: f32) bool {
    return gjk(a_pos, Collider{ .sphere = .{ .radius = a_rad } }, b_pos, Collider{ .sphere = .{ .radius = b_rad } });
}

pub fn gjk(a_pos: vec3, collider_a: Collider, b_pos: vec3, collider_b: Collider) bool {
    return gjk_inner(a_pos, collider_a, b_pos, collider_b, 24).collision;
}

pub fn gjk_inner(a_pos: vec3, collider_a: Collider, b_pos: vec3, collider_b: Collider, iters: usize) struct { collision: bool, simplex: Simplex } {
    // ref: https://winter.dev/articles/gjk-algorithm
    // ref: https://gist.github.com/vurtun/29727217c269a2fbf4c0ed9a1d11cb40
    var v = vec3.init(1, 0, 0);
    var a = collider_a.support(a_pos, v).sub(collider_b.support(b_pos, v.mul(-1)));
    var simplex = Simplex.new(&.{a});
    var d = a.mul(-1).normalize();

    var i: usize = 0;

    for (0..iters) |_| {
        i += 1;
        a = collider_a.support(a_pos, d).sub(collider_b.support(b_pos, d.mul(-1)));

        if (a.dot(d) < 0) {
            return .{ .collision = false, .simplex = simplex };
        }
        simplex.add(a);
        const contains_origin = next_simplex(&simplex, &d);
        d = d.normalize();
        if (contains_origin) {
            return .{ .collision = true, .simplex = simplex };
        }
    }

    std.debug.print("shouldn't get here.\n", .{});
    return .{ .collision = false, .simplex = simplex };
}

// TODO: https://winter.dev/articles/epa-algorithm

pub const ColliderType = enum {
    sphere,
    mesh,
};

pub const Collider = union(ColliderType) {
    sphere: struct { radius: f32 },
    mesh: struct { verts: []const vec3 },

    const Self = @This();
    fn support(self: Self, pos: vec3, dir: vec3) vec3 {
        const point = switch (self) {
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
        return pos.add(point);
    }
};

pub const Simplex = struct {
    points: BoundedArray(vec3, 4),

    const Self = @This();

    fn new(vs: []const vec3) Self {
        return Self{ .points = BoundedArray(vec3, 4).fromSlice(vs) catch unreachable };
    }

    fn add(self: *Self, v: vec3) void {
        self.points.insert(0, v) catch @panic("too many points");
    }
};

fn next_simplex(simplex: *Simplex, dir: *vec3) bool {
    switch (simplex.points.len) {
        2 => {
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
