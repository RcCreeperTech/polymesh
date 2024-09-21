const std = @import("std");
const rl = @import("raylib");

const EPSILON = std.math.pow(2, -52);
const EDGE_STACK = [512]u32;

const Allocator = std.mem.Allocator;

pub fn orient2d(ax: f32, ay: f32, bx: f32, by: f32, cx: f32, cy: f32) f32 {
    return (ay - cy) * (bx - cx) - (ax - cx) * (by - cy);
}

const Infinity = std.math.inf(f32);

pub const Delaunator = struct {
    // arrays that will store the triangulation graph
    coords: []f32,
    triangles: []u32,
    halfedges: []i32,

    // temporary arrays for tracking the edges of the advancing convex hull
    hullPrev: []u32,
    hullNext: []u32,
    hullTri: []u32,
    hullHash: []i32,
    hullStart: u32,

    // temporary arrays for sorting points
    ids: []u32,
    dists: []f64,

    cx: f32,
    cy: f32,

    pub fn init(allocator: std.mem.Allocator, points: []rl.Vector2) @This() {
        var self: @This() = undefined;

        self.coords = allocator.alloc(f32, points.len * 2);
        for (0..points.len) |i| {
            const p = points[i];
            self.coords[2 * i] = p.x;
            self.coords[2 * i + 1] = p.y;
        }

        // arrays that will store the triangulation graph
        const maxTriangles = @max(2 * points.len - 5, 0);
        self.triangles = allocator.alloc(u32, maxTriangles * 3);
        self.halfedges = allocator.alloc(i32, maxTriangles * 3);
        // temporary arrays for tracking the edges of the advancing convex hull
        self.hullPrev = allocator.alloc(u32, points.len); // edge to prev edge
        self.hullNext = allocator.alloc(u32, points.len); // edge to next edge
        self.hullTri = allocator.alloc(u32, points.len); // edge to adjacent triangle
        self.hullHash = allocator.alloc(i32, @ceil(@sqrt(points.len))); // angular edge hash

        // temporary arrays for sorting points
        self.ids = allocator.alloc(u32, points.len);
        self.dists = allocator.alloc(f64, points.len);

        self.update();

        return self;
    }

    pub fn update(self: *@This(), allocator: Allocator) void {
        var n = self.coords.length >> 1;

        // populate an array of point indices; calculate input data bbox
        const minX = Infinity;
        const minY = Infinity;
        const maxX = -Infinity;
        const maxY = -Infinity;

        for (0..n) |i| {
            const x = self.coords[2 * i];
            const y = self.coords[2 * i + 1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
            self.ids[i] = i;
        }
        const cx = (minX + maxX) / 2;
        const cy = (minY + maxY) / 2;

        var i_0: usize = undefined;
        var i_1: usize = undefined;
        var i_2: usize = undefined;

        // pick a seed point close to the center
        var minDist = Infinity;
        for (0..n) |i| {
            const d = dist(cx, cy, self.coords[2 * i], self.coords[2 * i + 1]);
            if (d < minDist) {
                i_0 = i;
                minDist = d;
            }
        }
        const i_0x = self.coords[2 * i_0];
        const i_0y = self.coords[2 * i_0 + 1];

        // find the point closest to the seed
        minDist = Infinity;
        for (0..n) |i| {
            if (i == i_0) continue;
            const d = dist(i_0x, i_0y, self.coords[2 * i], self.coords[2 * i + 1]);
            if (d < minDist and d > 0) {
                i_1 = i;
                minDist = d;
            }
        }
        var i_1x = self.coords[2 * i_1];
        var i_1y = self.coords[2 * i_1 + 1];

        var minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (0..n) |i| {
            if (i == i_0 or i == i_1) continue;
            const r = circumradius(
                i_0x,
                i_0y,
                i_1x,
                i_1y,
                self.coords[2 * i],
                self.coords[2 * i + 1],
            );
            if (r < minRadius) {
                i_2 = i;
                minRadius = r;
            }
        }
        var i_2x = self.coords[2 * i_2];
        var i_2y = self.coords[2 * i_2 + 1];

        if (minRadius == Infinity) {
            // order collinear points by dx (or dy if all x are identical)
            // and return the list as a hull
            for (0..n) |i| {
                self.dists[i] = (self.coords[2 * i] - self.coords[0]) || (self.coords[2 * i + 1] - self.coords[1]);
            }
            quicksort(self.ids, self.dists, 0, n - 1);
            const hull = allocator.alloc(u32, n);
            var j = 0;
            var d0 = -Infinity;
            for (0..n) |i| {
                const id = self.ids[i];
                const d = self.dists[id];
                if (d > d0) {
                    hull[j] = id;
                    j += 1;
                    d0 = d;
                }
            }
            self.hull = hull.subarray(0, j);
            // this.triangles = new Uint32Array(0);
            // this.halfedges = new Uint32Array(0);
            return;
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (orient2d(i_0x, i_0y, i_1x, i_1y, i_2x, i_2y) < 0) {
            const i = i_1;
            const x = i_1x;
            const y = i_1y;
            i_1 = i_2;
            i_1x = i_2x;
            i_1y = i_2y;
            i_2 = i;
            i_2x = x;
            i_2y = y;
        }

        const center = circumcenter(i_0x, i_0y, i_1x, i_1y, i_2x, i_2y);
        self.cx = center.x;
        self.cy = center.y;

        for (0..n) |i| {
            self.dists[i] = dist(self.coords[2 * i], self.coords[2 * i + 1], center.x, center.y);
        }

        // sort the points by distance from the seed triangle circumcenter
        quicksort(self.ids, self.dists, 0, n - 1);

        // set up the seed triangle as the starting hull
        self.hullStart = i_0;
        var hullSize = 3;

        self.hullNext[i_0] = i_1;
        self.hullNext[i_1] = i_2;
        self.hullNext[i_2] = i_0;

        self.hullPrev[i_2] = i_1;
        self.hullPrev[i_0] = i_2;
        self.hullPrev[i_1] = i_0;

        self.hullTri[i_0] = 0;
        self.hullTri[i_1] = 1;
        self.hullTri[i_2] = 2;

        self.hullHash.fill(-1);
        self.hullHash[self.hashKey(i_0x, i_0y)] = i_0;
        self.hullHash[self.hashKey(i_1x, i_1y)] = i_1;
        self.hullHash[self.hashKey(i_2x, i_2y)] = i_2;

        self.trianglesLen = 0;
        self.addTriangle(i_0, i_1, i_2, -1, -1, -1);

        var xp: f32 = undefined;
        var yp: f32 = undefined;
        for (0..self.ids.len) |k| {
            const i = self.ids[k];
            const x = self.coords[2 * i];
            const y = self.coords[2 * i + 1];

            // skip near-duplicate points
            if (k > 0 and @abs(x - xp) <= EPSILON and @abs(y - yp) <= EPSILON) continue;
            xp = x;
            yp = y;

            // skip seed triangle points
            if (i == i_0 or i == i_1 or i == i_2) continue;

            // find a visible edge on the convex hull using edge hash
            var start = 0;
            const key = self.hashKey(x, y);
            for (0..self.hullHash.len) |j| {
                start = self.hullHash[(key + j) % self.hullHash.len];
                if (start != -1 and start != self.hullNext[start]) break;
            }

            start = self.hullPrev[start];
            var e = start;
            var q = self.hullNext[e];
            while (orient2d(x, y, self.coords[2 * e], self.coords[2 * e + 1], self.coords[2 * q], self.coords[2 * q + 1]) >= 0) : (q = self.hullNext[e]) {
                e = q;
                if (e == start) {
                    e = -1;
                    break;
                }
            }
            if (e == -1) continue; // likely a near-duplicate point; skip it

            // add the first triangle from the point
            var t = self.addTriangle(e, i, self.hullNext[e], -1, -1, self.hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            self.hullTri[i] = self.legalize(t + 2);
            self.hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize += 1;

            // walk forward through the hull, adding more triangles and flipping recursively
            n = self.hullNext[e];
            q = self.hullNext[n];
            while (orient2d(x, y, self.coords[2 * n], self.coords[2 * n + 1], self.coords[2 * q], self.coords[2 * q + 1]) < 0) : (q = self.hullNext[n]) {
                t = self.addTriangle(n, i, q, self.hullTri[i], -1, self.hullTri[n]);
                self.hullTri[i] = self.legalize(t + 2);
                self.hullNext[n] = n; // mark as removed
                hullSize -= 1;
                n = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e == start) {
                q = self.hullNext[n];
                while (orient2d(x, y, self.coords[2 * q], self.coords[2 * q + 1], self.coords[2 * e], self.coords[2 * e + 1]) < 0) : (q = self.hullPrev[e]) {
                    t = self.addTriangle(q, i, e, -1, self.hullTri[e], self.hullTri[q]);
                    self.legalize(t + 2);
                    self.hullTri[q] = t;
                    self.hullNext[e] = e; // mark as removed
                    hullSize -= 1;
                    e = q;
                }
            }

            // update the hull indices
            self.hullStart = e;
            self.hullPrev[i] = e;
            self.hullNext[e] = i;
            self.hullPrev[n] = i;
            self.hullNext[i] = n;

            // save the two new edges in the hash table
            self.hullHash[self.hashKey(x, y)] = i;
            self.hullHash[self.hashKey(self.coords[2 * e], self.coords[2 * e + 1])] = e;
        }

        var e = self.hullStart;
        self.hull = allocator.alloc(u32, hullSize);
        for (0..self.hullStart) |i| {
            self.hull[i] = e;
            e = self.hullNext[e];
        }

        // trim typed triangle mesh arrays
        self.triangles = self.triangles.subarray(0, self.trianglesLen);
        self.halfedges = self.halfedges.subarray(0, self.trianglesLen);
    }

    fn hashKey(self: *@This(), x: f32, y: f32) usize {
        return @floor(pseudoAngle(x - self.cx, y - self.cy) * self.hashSize) % self.hashSize;
    }

    fn legalize(self: *@This(), a: usize) void {
        var i = 0;
        var ar = 0;

        // recursion eliminated with a fixed-size stack
        while (true) {
            const b = self.halfedges[a];

            // * if the pair of triangles doesn't satisfy the Delaunay condition
            // * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
            // * then do the same check/flip recursively for the new pair of triangles
            // *
            // *           pl                    pl
            // *          /||\                  /  \
            // *       al/ || \bl            al/    \a
            // *        /  ||  \              /      \
            // *       /  a||b  \    flip    /___ar___\
            // *     p0\   ||   /p1   =>   p0\---bl---/p1
            // *        \  ||  /              \      /
            // *       ar\ || /br             b\    /br
            // *          \||/                  \  /
            // *           pr                    pr
            // *
            const a0 = a - a % 3;
            ar = a0 + (a + 2) % 3;

            if (b == -1) { // convex hull edge
                if (i == 0) break;
                a = EDGE_STACK[--i];
                continue;
            }

            const b0 = b - b % 3;
            const al = a0 + (a + 1) % 3;
            const bl = b0 + (b + 2) % 3;

            const p0 = self.triangles[ar];
            const pr = self.triangles[a];
            const pl = self.triangles[al];
            const p1 = self.triangles[bl];

            const illegal = inCircle(
                self.coords[2 * p0],
                self.coords[2 * p0 + 1],
                self.coords[2 * pr],
                self.coords[2 * pr + 1],
                self.coords[2 * pl],
                self.coords[2 * pl + 1],
                self.coords[2 * p1],
                self.coords[2 * p1 + 1],
            );

            if (illegal) {
                self.triangles[a] = p1;
                self.triangles[b] = p0;

                const hbl = self.halfedges[bl];

                // edge swapped on the other side of the hull (rare); fix the halfedge reference
                if (hbl == -1) {
                    var e = self.hullStart;
                    if (self.hullTri[e] == bl) {
                        self.hullTri[e] = a;
                        break;
                    }
                    e = self.hullPrev[e];
                    while (e != self.hullStart) {
                        if (self.hullTri[e] == bl) {
                            self.hullTri[e] = a;
                            break;
                        }
                        e = self.hullPrev[e];
                    }
                }
                self.link(a, hbl);
                self.link(b, self.halfedges[ar]);
                self.link(ar, bl);

                const br = b0 + (b + 1) % 3;

                // don't worry about hitting the cap: it can only happen on extremely degenerate input
                if (i < EDGE_STACK.length) {
                    EDGE_STACK[i] = br;
                    i += 1;
                }
            } else {
                if (i == 0) break;
                i -= 1;
                a = EDGE_STACK[i];
            }
        }

        return ar;
    }

    fn link(self: *@This(), a: usize, b: usize) void {
        self.halfedges[a] = b;
        if (b != -1) self.halfedges[b] = a;
    }

    // add a new triangle given vertex indices and adjacent half-edge ids
    fn addTriangle(self: *@This(), i_0: usize, i_1: usize, i_2: usize, a: f32, b: f32, c: f32) usize {
        const t = self.trianglesLen;

        self.triangles[t] = i_0;
        self.triangles[t + 1] = i_1;
        self.triangles[t + 2] = i_2;

        self.link(t, a);
        self.link(t + 1, b);
        self.link(t + 2, c);

        self.trianglesLen += 3;

        return t;
    }
};

// monotonically increases with real angle, but doesn't need expensive trigonometry
fn pseudoAngle(dx: f32, dy: f32) f32 {
    const p = dx / (@abs(dx) + @abs(dy));
    return if (dy > 0) (3 - p) / 4 else (1 + p) / 4; // [0..1]
}

fn dist(ax: f32, ay: f32, bx: f32, by: f32) f32 {
    const dx = ax - bx;
    const dy = ay - by;
    return dx * dx + dy * dy;
}

fn inCircle(ax: f32, ay: f32, bx: f32, by: f32, cx: f32, cy: f32, px: f32, py: f32) f32 {
    const dx = ax - px;
    const dy = ay - py;
    const ex = bx - px;
    const ey = by - py;
    const fx = cx - px;
    const fy = cy - py;

    const ap = dx * dx + dy * dy;
    const bp = ex * ex + ey * ey;
    const cp = fx * fx + fy * fy;

    return dx * (ey * cp - bp * fy) -
        dy * (ex * cp - bp * fx) +
        ap * (ex * fy - ey * fx) < 0;
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

fn circumcenter(ax: f32, ay: f32, bx: f32, by: f32, cx: f32, cy: f32) rl.Vector2 {
    const dx = bx - ax;
    const dy = by - ay;
    const ex = cx - ax;
    const ey = cy - ay;

    const bl = dx * dx + dy * dy;
    const cl = ex * ex + ey * ey;
    const d = 0.5 / (dx * ey - dy * ex);

    const x = ax + (ey * bl - dy * cl) * d;
    const y = ay + (dx * cl - ex * bl) * d;

    return .{ .x = x, .y = y };
}

fn quicksort(ids: []u32, dists: []f64, left: usize, right: usize) void {
    _ = ids;
    _ = dists;
    _ = left;
    _ = right;
    // if (right - left <= 20) {
    //     for (let i = left + 1; i <= right; i++) {
    //         const temp = ids[i];
    //         const tempDist = dists[temp];
    //         let j = i - 1;
    //         while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
    //         ids[j + 1] = temp;
    //     }
    // } else {
    //     const median = (left + right) >> 1;
    //     let i = left + 1;
    //     let j = right;
    //     swap(ids, median, i);
    //     if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
    //     if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
    //     if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);
    //
    //     const temp = ids[i];
    //     const tempDist = dists[temp];
    //     while (true) {
    //         do i++; while (dists[ids[i]] < tempDist);
    //         do j--; while (dists[ids[j]] > tempDist);
    //         if (j < i) break;
    //         swap(ids, i, j);
    //     }
    //     ids[left + 1] = ids[j];
    //     ids[j] = temp;
    //
    //     if (right - i + 1 >= j - left) {
    //         quicksort(ids, dists, i, right);
    //         quicksort(ids, dists, left, j - 1);
    //     } else {
    //         quicksort(ids, dists, left, j - 1);
    //         quicksort(ids, dists, i, right);
    //     }
    // }
}

// fn swap(arr, i, j) {
//     const tmp = arr[i];
//     arr[i] = arr[j];
//     arr[j] = tmp;
// }

// function defaultGetX(p) {
//     return p[0];
// }
// function defaultGetY(p) {
//     return p[1];
// }
