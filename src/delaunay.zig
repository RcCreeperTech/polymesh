const std = @import("std");
const poly = @import("polymesh.zig");
const rl = @import("raylib");

const Vec2 = rl.Vector2;

const EPSILON = std.math.pow(2, -52);
const EDGE_STACK = [512]u32;
const Infinity = std.math.inf(f32);

const Allocator = std.mem.Allocator;

pub fn Delunator(comptime Mesh: type) type {
    // TODO: Type checking for a PolyMesh
    // TODO: Assert that vertex has a .pos field or is a @Vector(2, f32)...
    return struct {
        const Self = @This();

        // temporary arrays for tracking the edges of the advancing convex hull
        const ConvexHull = struct {
            start: Mesh.Vertex.Handle,
            // NOTE: Multiarraylist???
            prev: []Mesh.HalfEdge.Handle,
            next: []Mesh.HalfEdge.Handle,
            tri: []Mesh.Face.Handle,
            // Hash
            hash: std.HashMapUnmanaged(
                Mesh.Vertex.Handle,
                Mesh.HalfEdge.Handle,
                void, // FIXME: Context with pseudoAngle here
                std.hash_map.default_max_load_percentage,
            ),
        };

        mesh: *Mesh,
        hull: ConvexHull,
        // temporary arrays for sorting points
        ids: []Mesh.Vertex.Handle,
        dists: []f64,

        center: Vec2,

        pub fn init(allocator: Allocator, mesh: *Mesh) !Self {
            const num_points = mesh.verts.count();
            return .{
                .mesh = mesh,
                .hull = .{
                    .start = undefined,
                    .prev = try allocator.alloc(Mesh.HalfEdge.Handle, num_points),
                    .next = try allocator.alloc(Mesh.HalfEdge.Handle, num_points),
                    .tri = try allocator.alloc(Mesh.Face.Handle, num_points),
                    // angular edge hash
                    .hash = try allocator.alloc(Mesh.HalfEdge.Handle, @ceil(@sqrt(num_points))),
                },
                .ids = try allocator.alloc(Mesh.HalfEdge.Handle, num_points),
                .dists = try allocator.alloc(f64, num_points),
                .cy = undefined,
                .cx = undefined,
            };
        }

        pub fn triangulate(self: *Self) void {
            // TODO: Do nothing if there are not enough points
            var it = self.mesh.verts.iterator();
            // populate an array of point indices; calculate input data bbox
            self.center = compute_center: {
                var min = Vec2.init(-Infinity, -Infinity);
                var max = Vec2.init(Infinity, Infinity);
                var idx = 0;
                while (it.next()) |i| {
                    const vert: Mesh.Vertex.Handle = .from(i);
                    self.ids[idx] = vert;
                    idx += 1;
                    const pos: Vec2 = vert.deref(&self.mesh.verts).data.pos;
                    min = Vec2.min(min, pos);
                    max = Vec2.max(max, pos);
                }
                break :compute_center min.add(max).scale(0.5);
            };

            // Choose the 3 vertices that will comprise the initial triangle
            it.reset();
            var seed: [3]struct { handle: Mesh.Vertex.Handle, pos: Vec2 } = undefined;
            var min_val = Infinity;
            inline for (0..3) |idx| {
                while (it.next()) |handle| {
                    const vert: Mesh.Vertex.Handle = .from(handle);
                    const pos: Vec2 = vert.deref(&self.mesh.verts).data.pos;

                    const val = switch (idx) {
                        0 => pos.distance(self.center), // Find the point closest to the center
                        1 => pos.distance(seed[0].pos), // Find the point closest to seed[0]
                        2 => circumradius( // Find the point with min circumradius
                            seed[0].pos.x,
                            seed[0].pos.y,
                            seed[1].pos.x,
                            seed[1].pos.y,
                            pos.x,
                            pos.y,
                        ),
                        else => unreachable,
                    };
                    if (val < min_val) {
                        min_val = val;
                        seed[idx].handle = vert;
                    }
                }
                seed[idx].pos = seed[idx].handle.deref(&self.mesh.verts).data.pos;
            }
        }
        // if (minRadius == Infinity) {
        //            // order collinear points by dx (or dy if all x are identical)
        //            // and return the list as a hull
        //            for (0..n) |i| {
        //                self.dists[i] = (self.coords[2 * i] - self.coords[0]) || (self.coords[2 * i + 1] - self.coords[1]);
        //            }
        //            quicksort(self.ids, self.dists, 0, n - 1);
        //            const hull = allocator.alloc(u32, n);
        //            var j = 0;
        //            var d0 = -Infinity;
        //            for (0..n) |i| {
        //                const id = self.ids[i];
        //                const d = self.dists[id];
        //                if (d > d0) {
        //                    hull[j] = id;
        //                    j += 1;
        //                    d0 = d;
        //                }
        //            }
        //            self.hull = hull.subarray(0, j);
        //            // this.triangles = new Uint32Array(0);
        //            // this.halfedges = new Uint32Array(0);
        //            return;
        //        }
        //
        //        // swap the order of the seed points for counter-clockwise orientation
        //        if (orient2d(i_0x, i_0y, i_1x, i_1y, i_2x, i_2y) < 0) {
        //            const i = i_1;
        //            const x = i_1x;
        //            const y = i_1y;
        //            i_1 = i_2;
        //            i_1x = i_2x;
        //            i_1y = i_2y;
        //            i_2 = i;
        //            i_2x = x;
        //            i_2y = y;
        //        }

    };
}

fn orient2d(ax: f32, ay: f32, bx: f32, by: f32, cx: f32, cy: f32) f32 {
    return (ay - cy) * (bx - cx) - (ax - cx) * (by - cy);
}

fn circumradius(ax: f32, ay: f32, bx: f32, by: f32, cx: f32, cy: f32) f32 {
    const dx = bx - ax;
    const dy = by - ay;
    const ex = cx - ax;
    const ey = cy - ay;

    const bl = dx * dx + dy * dy;
    const cl = ex * ex + ey * ey;
    const d = 0.5 / (dx * ey - dy * ex);

    const x = (ey * bl - dy * cl) * d;
    const y = (dx * cl - ex * bl) * d;

    return x * x + y * y;
}
