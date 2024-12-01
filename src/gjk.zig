const rl = @import("raylib");
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

    pub fn mag(self: Self) f32 {
        return @sqrt(self.dot(self));
    }

    pub fn normalize(self: Self) Self {
        const len = self.mag();
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

const epa_max_points = 32;
const EpaVec3Arr = std.BoundedArray(vec3, epa_max_points);
fn getFaceNormals(polytope: EpaVec3Arr, faces: std.BoundedArray(usize, epa_max_points * 3)) struct { normals: EpaVec3Arr, distances: [epa_max_points]f32, min_face: usize } {
    var min_triangle: usize = 0;
    var min_dist = std.math.floatMax(f32);
    var normals = EpaVec3Arr.init(0) catch unreachable;
    var distances = [_]f32{0} ** epa_max_points;

    var i: usize = 0;
    while (i < faces.len) : (i += 3) {
        const a = polytope.get(faces.get(i));
        const b = polytope.get(faces.get(i + 1));
        const c = polytope.get(faces.get(i + 2));

        var normal = b.sub(a).cross(c.sub(a)).normalize();
        var distance = normal.dot(a);

        if (distance < 0) {
            normal = normal.mul(-1);
            distance = -distance;
        }

        normals.append(normal) catch @panic("normals");
        distances[normals.len - 1] = distance;
        if (distance < min_dist) {
            min_triangle = i / 3;
            min_dist = distance;
        }
    }

    return .{
        .normals = normals,
        .distances = distances,
        .min_face = min_triangle,
    };
}

fn addIfUniqueEdge(
    edges: *std.BoundedArray(std.meta.Tuple(&.{ usize, usize }), epa_max_points),
    faces: std.BoundedArray(usize, epa_max_points * 3),
    a: usize,
    b: usize,
) void {
    const reverse = lbl: {
        for (edges.slice(), 0..) |e, i| {
            if (std.meta.eql(e, .{ faces.get(b), faces.get(a) })) {
                break :lbl i;
            }
        }
        break :lbl null;
    };

    if (reverse) |r| {
        _ = edges.orderedRemove(r);
    } else {
        edges.append(.{ faces.get(a), faces.get(b) }) catch @panic("edges");
    }
}

const Polytope = struct {
    verts: []const vec3,
    faces: [][3]usize,
};

const Dimension = enum { x, y, z };
fn ltDim(dim: Dimension, lhs: vec3, rhs: vec3) bool {
    return switch (dim) {
        .x => lhs.x < rhs.x,
        .y => lhs.y < rhs.y,
        .z => lhs.z < rhs.z,
    };
}

const Line = struct { bc_dir: vec3, b: vec3 };
fn distToLine(l: Line, a: vec3) f32 {
    // https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
    const v = a.sub(l.b);
    const t = v.dot(l.bc_dir);
    return l.b.add(l.bc_dir.mul(t)).sub(a).mag();
}

/// argmax using a `key` mapping function rather than `lessThan`
pub fn argMax2(
    comptime T: type,
    items: []const T,
    context: anytype,
    comptime key: fn (@TypeOf(context), val: T) f32,
) ?struct { f32, usize } {
    if (items.len == 0) {
        return null;
    }

    var largest = key(context, items[0]);
    var largest_index: usize = 0;
    for (items[1..], 1..) |item, i| {
        const k = key(context, item);
        if (k > largest) {
            largest = k;
            largest_index = i;
        }
    }

    return .{ largest, largest_index };
}

fn absdot(a: vec3, b: vec3) f32 {
    return @abs(a.dot(b));
}

const Plane = struct { point: vec3, normal: vec3 };
fn planeDist(plane: Plane, point: vec3) f32 {
    return point.sub(plane.point).dot(plane.normal);
}

fn absPlaneDist(plane: Plane, point: vec3) f32 {
    return @abs(planeDist(plane, point));
}

fn triNormal(triangle: [3]vec3) vec3 {
    const p = triangle[0];
    const q = triangle[1];
    const r = triangle[2];
    return r.sub(p).cross(q.sub(p)).normalize();
}

const Face = [3]usize;

const TriMesh = struct {
    const HalfEdge = struct {
        vertex: usize,
        face: usize,
        next: [2]usize,
        opposite: ?[2]usize,
    };

    allocator: std.mem.Allocator,
    vertices: []const vec3,
    faces: std.ArrayList(Face),

    // map edge (vertex indices) to HalfeEdge
    hemap: std.AutoHashMap([2]usize, HalfEdge),

    const Self = @This();

    fn deinit(self: *Self) void {
        self.hemap.deinit();
        self.faces.deinit();
        //self.allocator.free(self.vertices);
    }

    fn addFace(self: *Self, verts: [3]usize) !void {
        std.debug.assert(new_face: {
            for (self.faces.items, 0..) |face, i| {
                if (std.meta.eql(face, verts)) {
                    std.debug.print("already have {any}={any} at {}\n", .{ verts, face, i });
                    break :new_face false;
                }
            }
            break :new_face true;
        });
        const face_id = self.faces.items.len;

        try self.faces.append(verts);

        for (0..3) |j| {
            const u = verts[j];
            const v = verts[(j + 1) % 3];
            const w = verts[(j + 2) % 3];

            try self.hemap.putNoClobber(.{ u, v }, HalfEdge{
                .vertex = verts[j],
                .face = face_id,
                .next = .{ v, w },
                .opposite = opp: {
                    if (self.hemap.getPtr(.{ v, u })) |opp| {
                        opp.opposite = .{ u, v };
                        break :opp .{ v, u };
                    } else {
                        break :opp null;
                    }
                },
            });
        }
    }

    fn build(allocator: std.mem.Allocator, vertices: []const vec3, faces: []const [3]usize) !Self {
        const hemap = std.AutoHashMap([2]usize, HalfEdge).init(allocator);
        const hefaces = try std.ArrayList(Face).initCapacity(allocator, faces.len);

        var r = Self{
            .allocator = allocator,
            .vertices = vertices,
            .faces = hefaces,
            .hemap = hemap,
        };

        for (faces) |face| {
            try r.addFace(face);
        }
        return r;
    }
};

// fn computeHalfEdges(faces: [][3]usize) void {
//     std.debug.assert(verts: {
//         var m = 0;
//         for (faces) |face| {
//             for (face) |i| {
//                 m = @max(m, i);
//             }
//         }
//         break :verts m + 1;
//     } == out.len);
//
//     for (faces, 0..) |face, face_i| {
//         for (.{ .{ face[0], face[1] }, .{ face[1], face[2] }, .{ face[2], face[0] } }) |edge| {
//             const u, const v = edge;
//             const next_i = face[(face_vert - 1) % 3];
//             for (faces, 0..) |f, j| {
//                 if ((f[0] == next_i and f[1] == i) or
//                     (f[1] == next_i and f[2] == i) or
//                     (f[2] == next_i and f[0] == i))
//                 {}
//             }
//             half_edges[i] = HalfEdge{
//                 .vertex = i,
//                 .face = face_i,
//                 .next = next_i,
//                 .opposite = null,
//             };
//         }
//     }
// }

fn faceEdges(face: Face) [3][2]usize {
    return .{
        .{ face[0], face[1] },
        .{ face[1], face[2] },
        .{ face[2], face[0] },
    };
}

fn reorderEdgeLoop(edges: [][2]usize) void {
    for (0..edges.len) |i| {
        for (i..edges.len) |j| {
            if (edges[i][1] == edges[j][0] and j != i + 1) {
                std.mem.swap([2]usize, &edges[i + 1], &edges[j]);
                break;
            }
        }
        std.debug.assert(edges[i][1] == edges[(i + 1) % edges.len][0]);
    }
}

pub fn quickhull(allocator: std.mem.Allocator, points: []const vec3) !?Polytope {
    if (points.len <= 3) {
        return null;
    }
    var outside_points = try std.DynamicBitSet.initFull(allocator, points.len);
    defer outside_points.deinit();

    var polytope = Polytope{ .verts = &.{}, .faces = &.{} };

    var mesh: TriMesh = undefined;
    // (1) Find maximum simplex
    {
        var simplex: [4]usize = undefined;
        simplex_line: {
            inline for (std.meta.fields(Dimension)) |d| {
                const mini = std.sort.argMin(vec3, points, @field(Dimension, d.name), ltDim).?;
                const maxi = std.sort.argMax(vec3, points, @field(Dimension, d.name), ltDim).?;
                if (mini != maxi) {
                    simplex[0] = mini;
                    simplex[1] = maxi;

                    break :simplex_line;
                }
            }
            // degenerate case, min = max along all dimensions, 0d
            return null;
        }

        { // find point farthest from line
            const p = points[simplex[0]];
            const q = points[simplex[1]];
            const pq_dir = p.sub(q).normalize();
            const linedist, const i = argMax2(vec3, points, Line{ .bc_dir = pq_dir, .b = p }, distToLine).?;
            if (linedist == 0) {
                // degenerate 1d case
                return null;
            }

            simplex[2] = i;
        }

        { // find point farthest from plane
            const p = points[simplex[0]];
            const q = points[simplex[1]];
            const r = points[simplex[2]];
            const planedist, const i = argMax2(vec3, points, Plane{ .point = p, .normal = triNormal(.{ p, q, r }) }, absPlaneDist).?;

            if (planedist == 0) {
                // degenerate 2d case
                return null;
            }

            simplex[3] = i;
        }

        mesh = try TriMesh.build(allocator, try allocator.dupe(vec3, points), &.{
            .{ 0, 1, 2 },
            .{ 0, 3, 1 },
            .{ 0, 2, 3 },
            .{ 1, 3, 2 },
        });
    }
    defer mesh.deinit();

    // (2) remove inner points, find max outer points for each face
    var farthest_from_face = try allocator.alloc(std.meta.Tuple(&.{ f32, usize }), mesh.faces.items.len);
    defer allocator.free(farthest_from_face);
    @memset(farthest_from_face, .{ 0, 0 });

    {
        var iter = outside_points.iterator(.{});

        while (iter.next()) |i| {
            const point = points[i];
            var outside = false;
            for (0..mesh.faces.items.len) |j| {
                const face = mesh.faces.items[j];
                const p = points[face[0]];
                const q = points[face[1]];
                const r = points[face[2]];
                const dist = planeDist(.{ .point = p, .normal = triNormal(.{ p, q, r }) }, point);
                if (dist > 0) {
                    outside = true;
                }
                if (dist > farthest_from_face[j].@"0") {
                    farthest_from_face[j] = .{ dist, i };
                }
            }
            if (!outside) {
                rl.drawSphere(@bitCast(point), 0.1, rl.Color.gray);
                outside_points.unset(i);
            }
        }
    }

    std.debug.print("built mesh\n", .{});

    // (3) expand to farthest point for each face
    for (farthest_from_face, 0..) |f, fi| {
        const dist, const farthest_point_i = f;
        const farthest_point = points[farthest_point_i];
        if (dist > 0) {
            var edgeloop = std.ArrayList([2]usize).init(allocator);
            defer edgeloop.deinit();

            var horizon_facet = std.ArrayList(usize).init(allocator);
            try horizon_facet.append(fi);
            defer horizon_facet.deinit();

            @import("debug.zig").debugMesh(@ptrCast(mesh.vertices), mesh.faces.items);
            var fjj: usize = 0;
            while (fjj < horizon_facet.items.len) : (fjj += 1) {
                const fj = horizon_facet.items[fjj];
                for (faceEdges(mesh.faces.items[fj])) |edge| {
                    const a, const b = edge;
                    if (mesh.hemap.get(.{ b, a })) |opp_he| {
                        const already_visited = std.mem.indexOf(usize, horizon_facet.items, &.{opp_he.face}) != null;
                        if (!already_visited) {
                            const face = mesh.faces.items[opp_he.face];
                            const p1 = points[face[0]];
                            const p2 = points[face[1]];
                            const p3 = points[face[2]];
                            const norm = triNormal(.{ p1, p2, p3 });
                            if (norm.dot(farthest_point.sub(p1)) > 0) {
                                // face is visible from point
                                try horizon_facet.append(opp_he.face);
                            } else {
                                // not visible from point, so the edge we crossed is a horizon edge
                                try edgeloop.append(edge);
                            }
                        }
                    }
                }

                @import("debug.zig").highlightFace(fj);
            }

            reorderEdgeLoop(edgeloop.items);
            for (horizon_facet.items) |face| {
                for (faceEdges(mesh.faces.items[face])) |edge| {
                    _ = mesh.hemap.remove(edge);
                }
                mesh.faces.items[face] = .{ 0, 0, 0 };
                //maybe??
                //farthest_from_face[face].@"1" = 0;
            }

            for (edgeloop.items) |edge| {
                try mesh.addFace(.{ edge[0], edge[1], farthest_point_i });
            }

            if (@import("debug.zig").isDebugIter()) {
                std.debug.print("{any}\n", .{horizon_facet.items});
                rl.drawSphere(@bitCast(farthest_point), 0.4, rl.Color.red);
            }
        }
    }

    @import("debug.zig").debugMesh(@ptrCast(mesh.vertices), mesh.faces.items);

    // (4) return result
    polytope.verts = mesh.vertices;
    polytope.faces = try mesh.faces.toOwnedSlice();
    return polytope;
}

pub fn epa(simplex: Simplex, a_pos: vec3, collider_a: Collider, b_pos: vec3, collider_b: Collider) struct { normal: vec3, penetration_depth: f32 } {
    var polytope = EpaVec3Arr.fromSlice(simplex.points.slice()) catch unreachable;
    var faces = std.BoundedArray(usize, epa_max_points * 3).fromSlice(&.{
        0, 1, 2,
        0, 3, 1,
        0, 2, 3,
        1, 3, 2,
    }) catch unreachable;

    var n = getFaceNormals(polytope, faces);

    var min_norm = vec3.init(0, 0, 0);
    var min_dist = std.math.floatMax(f32);

    while (min_dist == std.math.floatMax(f32)) {
        min_norm = n.normals.get(n.min_face);
        min_dist = n.distances[n.min_face];

        const support = collider_a.support(a_pos, min_norm).sub(collider_b.support(b_pos, min_norm.mul(-1)));
        const s_dist = min_norm.dot(support);

        if (@abs(s_dist - min_dist) > 0.001) {
            min_dist = std.math.floatMax(f32);

            var unique_edges = std.BoundedArray(std.meta.Tuple(&.{ usize, usize }), epa_max_points).init(0) catch unreachable;
            {
                var i: usize = 0;
                while (i < n.normals.len) : (i += 1) {
                    if (n.normals.get(i).dot(support) > 0) {
                        const f = i * 3;

                        addIfUniqueEdge(&unique_edges, faces, f, f + 1);
                        addIfUniqueEdge(&unique_edges, faces, f + 1, f + 2);
                        addIfUniqueEdge(&unique_edges, faces, f + 2, f);

                        faces.set(f + 2, faces.pop());
                        faces.set(f + 1, faces.pop());
                        faces.set(f, faces.pop());

                        n.normals.set(i, n.normals.pop());

                        i -= 1;
                    }
                }
            }

            var new_faces = std.BoundedArray(usize, epa_max_points * 3).init(0) catch unreachable;

            for (unique_edges.slice()) |e| {
                new_faces.append(e.@"0") catch @panic("new face");
                new_faces.append(e.@"1") catch @panic("new face");
                new_faces.append(polytope.len) catch @panic("new face");
            }

            polytope.append(support) catch @panic("polytope");

            var new = getFaceNormals(polytope, new_faces);

            var old_min_distance = std.math.floatMax(f32);
            {
                var i: usize = 0;
                while (i < n.normals.len) : (i += 1) {
                    if (n.distances[i] < old_min_distance) {
                        old_min_distance = n.distances[i];
                        n.min_face = i;
                    }
                }
            }

            if (new.distances[new.min_face] < old_min_distance) {
                n.min_face = new.min_face + n.normals.len;
            }

            faces.appendSlice(new_faces.slice()) catch @panic("appending new faces");
            n.normals.appendSlice(new.normals.slice()) catch @panic("appending new normals");
        }
    }

    return .{
        .normal = min_norm,
        .penetration_depth = min_dist + 0.001,
    };
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
