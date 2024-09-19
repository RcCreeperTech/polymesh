const std = @import("std");
const rl = @import("raylib");

const EPSILON = std.math.pow(2, -52);
const EDGE_STACK = [512]u32;


const orient2d = @import("./robust_orient2d_translated.zig").orient2d;

const Infinity = std.math.inf(f32);

pub const Delaunator = struct {
    // arrays that will store the triangulation graph
    coords: []f32,
    triangles: []u32,
    halfedges: []i32,

    // temporary arrays for tracking the edges of the advancing convex hull
    hashSize: usize,
    hullPrev: []u32,
    hullNext: []u32,
    hullTri: []u32,
    hullHas: []i32,

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
        self.hashSize = @ceil(@sqrt(points.len));
        self.hullPrev = allocator.alloc(u32, points.len); // edge to prev edge
        self.hullNext = allocator.alloc(u32, points.len); // edge to next edge
        self.hullTri = allocator.alloc(u32, points.len); // edge to adjacent triangle
        self.hullHash = allocator.alloc(i32, self.hashSize); // angular edge hash

        // temporary arrays for sorting points
        self.ids = allocator.alloc(u32, points.len);
        self.dists = allocator.alloc(f64, points.len);

        self.update();

        return self;
    }

    pub fn update(self: *@This()) void {
        const n = self.coords.length >> 1;

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
                i0 = i;
                minDist = d;
            }
        }
        const i_0x = self.coords[2 * i0];
        const i_0y = self.coords[2 * i0 + 1];

        // find the point closest to the seed
        minDist = Infinity;
        for (0..n) |i| {
            if (i == i0) continue;
            const d = dist(i_0x, i_0y, self.coords[2 * i], self.coords[2 * i + 1]);
            if (d < minDist and d > 0) {
                i1 = i;
                minDist = d;
            }
        }
        var i_1x = self.coords[2 * i1];
        var i_1y = self.coords[2 * i1 + 1];

        var minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (0..n) |i| {
            if (i == i0 or i == i1) continue;
            const r = circumradius(i_0x, i_0y, i_1x, i_1y, self.coords[2 * i], self.coords[2 * i + 1]);
            if (r < minRadius) {
                i2 = i;
                minRadius = r;
            }
        }
        var i2x = self.coords[2 * i2];
        var i2y = self.coords[2 * i2 + 1];

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
            this.hull = hull.subarray(0, j);
            // this.triangles = new Uint32Array(0);
            // this.halfedges = new Uint32Array(0);
            return;
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0) {
            const i = i1;
            const x = i1x;
            const y = i1y;
            i1 = i2;
            i1x = i2x;
            i1y = i2y;
            i2 = i;
            i2x = x;
            i2y = y;
        }

        const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
        self.cx = center.x;
        self.cy = center.y;

        for (0..n) |i| {
            self.dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
        }

        // sort the points by distance from the seed triangle circumcenter
        quicksort(self.ids, self.dists, 0, n - 1);

        // set up the seed triangle as the starting hull
        self.hullStart = i0;
        var hullSize = 3;

        self.hullNext[i0] = i1;
        self.hullNext[i1] = i2;
        self.hullNext[i2] = i0;

        self.hullPrev[i2] = i1;
        self.hullPrev[i0] = i2;
        self.hullPrev[i1] = i0;

        hullTri[i0] = 0;
        hullTri[i1] = 1;
        hullTri[i2] = 2;

        hullHash.fill(-1);
        hullHash[self.hashKey(i0x, i0y)] = i0;
        hullHash[self.hashKey(i1x, i1y)] = i1;
        hullHash[self.hashKey(i2x, i2y)] = i2;

        this.trianglesLen = 0;
        self.addTriangle(i0, i1, i2, -1, -1, -1);

        var xp: f32 = undefined;
        var yp: f32 = undefined;
        for (0..self.ids.len) |k| {
            const i = self.ids[k];
            const x = coords[2 * i];
            const y = coords[2 * i + 1];

            // skip near-duplicate points
            if (k > 0 and @abs(x - xp) <= EPSILON and @abs(y - yp) <= EPSILON) continue;
            xp = x;
            yp = y;

            // skip seed triangle points
            if (i == i0 or i == i1 or i == i2) continue;

            // find a visible edge on the convex hull using edge hash
            var start = 0;
            var key = self.hashKey(x, y);
            for (0..self.hashSize) |j| {
                start = self.hullHash[(key + j) % self.hashSize];
                if (start != -1 and start != self.hullNext[start]) break;
            }

            start = self.hullPrev[start];
            var e = start;
            var q = hullNext[e];
            while (orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]) >= 0) : (q = hullNext[e]) {
                e = q;
                if (e == start) {
                    e = -1;
                    break;
                }
            }
            if (e == -1) continue; // likely a near-duplicate point; skip it

            // add the first triangle from the point
            var t = self.addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i] = self.legalize(t + 2);
            hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize += 1;

            // walk forward through the hull, adding more triangles and flipping recursively
            var n = hullNext[e];
            q = hullNext[n];
            while (orient2d(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1]) < 0) : (q = hullNext[n]) {
                t = self.addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                hullTri[i] = self.legalize(t + 2);
                hullNext[n] = n; // mark as removed
                hullSize -= 1;
                n = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e == start) {
                q = hullNext[n];
                while (orient2d(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]) < 0) : (q = hullPrev[e]) {
                    t = self.addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                    self.legalize(t + 2);
                    hullTri[q] = t;
                    hullNext[e] = e; // mark as removed
                    hullSize -= 1;
                    e = q;
                }
            }

            // update the hull indices
            self.hullStart = e;
            hullPrev[i] = e;
            hullNext[e] = i;
            hullPrev[n] = i;
            hullNext[i] = n;

            // save the two new edges in the hash table
            hullHash[self.hashKey(x, y)] = i;
            hullHash[self.hashKey(coords[2 * e], coords[2 * e + 1])] = e;
        }

        self.hull = allocator.alloc(u32, hullSize);
        for (0..self.hullStart) |i| {
            this.hull[i] = e;
            e = hullNext[e];
        }

        // trim typed triangle mesh arrays
        this.triangles = self.triangles.subarray(0, this.trianglesLen);
        this.halfedges = self.halfedges.subarray(0, this.trianglesLen);
    }

    fn hashKey(self: *@This(), x: f32, y: f32) usize {
        return @floor(pseudoAngle(x - self.cx, y - self.cy) * self.hashSize) % self.hashSize;
    }

    fn legalize(a) void {
        var i = 0;
        var ar = 0;

        // recursion eliminated with a fixed-size stack
        while (true) {
            const b = halfedges[a];

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

            const p0 = triangles[ar];
            const pr = triangles[a];
            const pl = triangles[al];
            const p1 = triangles[bl];

            const illegal = inCircle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl], coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]);

            if (illegal) {
                triangles[a] = p1;
                triangles[b] = p0;

                const hbl = halfedges[bl];

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
                self.link(b, halfedges[ar]);
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
    fn addTriangle(self: *@This(), i_0: usize, i_1: usize, i_2: usize, a: f32, b: f32, c: f32) void {
        _ = i_0; // autofix
        _ = i_1; // autofix
        _ = i_2; // autofix
        const t = self.trianglesLen;

        self.triangles[t] = i0;
        self.triangles[t + 1] = i1;
        self.triangles[t + 2] = i2;

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

fn quicksort(ids, dists, left, right) {
    if (right - left <= 20) {
        for (let i = left + 1; i <= right; i++) {
            const temp = ids[i];
            const tempDist = dists[temp];
            let j = i - 1;
            while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
            ids[j + 1] = temp;
        }
    } else {
        const median = (left + right) >> 1;
        let i = left + 1;
        let j = right;
        swap(ids, median, i);
        if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
        if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
        if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

        const temp = ids[i];
        const tempDist = dists[temp];
        while (true) {
            do i++; while (dists[ids[i]] < tempDist);
            do j--; while (dists[ids[j]] > tempDist);
            if (j < i) break;
            swap(ids, i, j);
        }
        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left) {
            quicksort(ids, dists, i, right);
            quicksort(ids, dists, left, j - 1);
        } else {
            quicksort(ids, dists, left, j - 1);
            quicksort(ids, dists, i, right);
        }
    }
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
