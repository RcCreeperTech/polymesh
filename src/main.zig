const std = @import("std");

const rl = @import("raylib");
const mesh = @import("mesh.zig");

const assert = std.debug.assert;

pub const std_options = .{
    // Set the log level to info
    .log_level = .debug,
};

const Line = struct {
    a: f32,
    b: f32,
    c: f32,

    pub fn evaluate(self: *@This(), x: f32, y: f32) f32 {
        return self.a * x + self.b * y + self.c;
    }
};

const Voronoi = struct {
    const log = std.log.scoped(.voronoi);

    const PolyMesh = mesh.GenericPolyMesh(rl.Vector2, void, void);

    const Beachline = struct {
        const Arc = struct {
            focus: rl.Vector2,
            event: ?*EventQueue.Event,
            boundary_left: ?*Boundary,
            boundary_right: ?*Boundary,

            pub const default: @This() = .{
                .focus = undefined,
                .event = null,
                .boundary_left = null,
                .boundary_right = null,
            };

            fn evaluate(self: *const @This(), x: f32, directrix_y: f32) f32 {
                const focus = self.focus;
                const diff_x = focus.x - x;
                const diff_y = focus.y - directrix_y;
                var result = std.math.pow(f32, diff_x, 2) / (2 * diff_y);
                result += (focus.y + directrix_y) / 2;
                return result;
            }

            fn compare(x: f32, item: @This()) std.math.Order {
                return std.math.order(item.focus.x, x);
            }

            pub fn format(
                self: @This(),
                comptime _: []const u8,
                _: std.fmt.FormatOptions,
                writer: anytype,
            ) !void {
                if (self.boundary_left != null)
                    try writer.print("->", .{})
                else
                    try writer.print("|-", .{});

                try writer.print("{{{d}, {d}}}", .{ self.focus.x, self.focus.y });
                if (self.event != null) try writer.print("!", .{});

                if (self.boundary_right != null)
                    try writer.print("<-", .{})
                else
                    try writer.print("-|", .{});
            }
        };
        const Boundary = struct {
            start: rl.Vector2,
            direction: rl.Vector2,
            mesh_hedge: *PolyMesh.HalfEdge,
        };

        arcs: std.ArrayListUnmanaged(Arc),
        boundary_pool: std.heap.MemoryPool(Boundary),

        pub fn init(allocator: std.mem.Allocator) @This() {
            return .{
                .arcs = .empty,
                .boundary_pool = .init(allocator),
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.arcs.deinit(allocator);
            self.boundary_pool.deinit();
        }

        pub fn findOverlappingArcInfo(self: *@This(), point: rl.Vector2) struct {
            usize, // arc_idx
            f32, // yval
        } {
            log.info("Searcing Point {{{d}, {d}}}", .{ point.x, point.y });
            // Boundary Conditions
            if (self.arcs.items.len == 0) return .{ 0, 0 };

            const index_of_first_arc_greater_than_or_equal_to_point = std.sort.lowerBound(Arc, self.arcs.items, point.x, Arc.compare);
            if (index_of_first_arc_greater_than_or_equal_to_point == 0) {
                log.info("Hit lower boundary", .{});
                const arc = self.getArc(0);
                const y = arc.evaluate(point.x, point.y);
                return .{ 0, y };
            } else if (index_of_first_arc_greater_than_or_equal_to_point == self.arcs.items.len) {
                log.info("Hit upper boundary", .{});
                const arc = self.arcs.getLast();
                const y = arc.evaluate(point.x, point.y);
                return .{ self.arcs.items.len - 1, y };
            } else {
                log.info("Hit middle boundary", .{});
                const lower_idx = index_of_first_arc_greater_than_or_equal_to_point - 1;
                const higher_idx = index_of_first_arc_greater_than_or_equal_to_point;

                const ylow = self.getArc(lower_idx).evaluate(point.x, point.y);
                const yhigh = self.getArc(higher_idx).evaluate(point.x, point.y);

                if (yhigh > ylow) {
                    return .{ higher_idx, yhigh };
                } else {
                    return .{ lower_idx, ylow };
                }
            }
        }

        pub fn getArc(self: *const @This(), index: usize) Arc {
            return self.arcs.items[index];
        }

        pub fn getArcPtr(self: *const @This(), index: usize) *Arc {
            return &self.arcs.items[index];
        }

        pub fn getArcPtrOrNull(self: *const @This(), index: usize) ?*Arc {
            return if (index >= self.arcs.items.len or index < 0) null else &self.arcs.items[index];
        }

        pub fn createBoundary(self: *@This(), boundary: Boundary) !*Boundary {
            const result = try self.boundary_pool.create();
            result.* = boundary;
            return result;
        }

        pub fn format(
            self: @This(),
            comptime _: []const u8,
            _: std.fmt.FormatOptions,
            writer: anytype,
        ) !void {
            try writer.print("[", .{});
            for (self.arcs.items, 0..) |arc, i| {
                try writer.print("{any}", .{arc});
                if (i != self.arcs.items.len - 1) try writer.print(",", .{});
            }
            return writer.print("]", .{});
        }
    };
    const EventQueue = struct {
        const Event = union(enum) {
            const PriorityQueue = std.PriorityQueue(*@This(), void, lessThan);
            site: rl.Vector2,
            circle: struct {
                cancelled: bool = false,
                center: rl.Vector2,
                radius: f32,
                arc: *Beachline.Arc,
            },

            pub fn sortHandle(self: *const @This()) f32 {
                return switch (self.*) {
                    .site => |pos| pos.y,
                    .circle => |circle| circle.center.y + circle.radius,
                };
            }

            pub fn lessThan(_: void, a: *@This(), b: *@This()) std.math.Order {
                return std.math.order(a.sortHandle(), b.sortHandle());
            }

            pub fn format(
                self: @This(),
                comptime _: []const u8,
                _: std.fmt.FormatOptions,
                writer: anytype,
            ) !void {
                return switch (self) {
                    .site => |s| try writer.print("site {{{d}, {d}}}", .{ s.x, s.y }),
                    .circle => |c| try writer.print("circle [{any}] center: {{{d}, {d}}} radius: {d}", .{
                        c.cancelled,
                        c.center.x,
                        c.center.y,
                        c.radius,
                    }),
                };
            }
        };

        queue: Event.PriorityQueue,
        event_pool: std.heap.MemoryPool(Event),

        pub fn init(allocator: std.mem.Allocator) @This() {
            return .{
                .queue = .init(allocator, {}),
                .event_pool = .init(allocator),
            };
        }

        pub fn deinit(self: *@This()) void {
            self.queue.deinit();
            self.event_pool.deinit();
        }

        pub fn count(self: *@This()) usize {
            return self.queue.count();
        }

        pub fn add(self: *@This(), elem: Event) !void {
            const e = try self.event_pool.create();
            e.* = elem;
            try self.queue.add(e);
        }

        pub fn addTracking(self: *@This(), elem: Event) !*Event {
            const e = try self.event_pool.create();
            e.* = elem;

            try self.queue.ensureUnusedCapacity(1);

            const start_index = self.queue.items.len - 1;
            self.queue.items.len += 1;
            self.queue.items[start_index] = e;

            const child = self.queue.items[start_index];
            var child_index = start_index;
            while (child_index > 0) {
                const parent_index = ((child_index - 1) >> 1);
                const parent = self.queue.items[parent_index];
                if (Event.lessThan(self.queue.context, child, parent) != .lt) break;
                self.queue.items[child_index] = parent;
                child_index = parent_index;
            }
            self.queue.items[child_index] = child;

            return self.queue.items[child_index];
        }

        pub fn removeOrNull(self: *@This()) ?*Event {
            return self.queue.removeOrNull();
        }
    };

    directrix: f32,
    mesh: PolyMesh,
    beachline: Beachline,
    event_queue: EventQueue,
    clip_rect: rl.Rectangle,

    pub fn init(allocator: std.mem.Allocator) @This() {
        return .{
            .directrix = 0,
            .mesh = .init(allocator),
            .beachline = .init(allocator),
            .event_queue = .init(allocator),
            .clip_rect = undefined, // FIXME:
        };
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        self.beachline.deinit(allocator);
        self.mesh.deinit();
        self.event_queue.deinit();
    }

    pub fn generate(self: *@This(), allocator: std.mem.Allocator, points: []rl.Vector2, optional_clip_rect: ?rl.Rectangle) !void {
        self.clip_rect = gen_clip_rect: {
            if (optional_clip_rect) |rect| {
                break :gen_clip_rect rect;
            } else {
                const aabb_padding = 20;
                const infinity = std.math.inf(f32);
                var min, var max = .{
                    rl.Vector2.init(infinity, infinity),
                    rl.Vector2.init(-infinity, -infinity),
                };

                for (points) |point| {
                    min = rl.Vector2.min(min, point);
                    max = rl.Vector2.max(max, point);
                }
                min = min.subtractValue(aabb_padding);
                max = max.addValue(aabb_padding);
                break :gen_clip_rect .{
                    .x = min.x,
                    .y = min.y,
                    .width = max.x - min.x,
                    .height = max.y - min.y,
                };
            }
        };

        for (points) |point| {
            if (rl.checkCollisionPointRec(point, self.clip_rect)) {
                try self.event_queue.add(
                    .{ .site = point },
                );
            }
        }
        // TODO: Resize the clip rect

        while (self.event_queue.count() > 0) try self.processEvent(allocator);

        // TODO: Terminate remaining edges on the bounding box
    }

    pub fn processEvent(self: *@This(), allocator: std.mem.Allocator) !void {
        if (self.event_queue.removeOrNull()) |event| {
            log.info("Processing Event {any} ", .{event});
            switch (event.*) {
                .site => |site| try self.siteEvent(allocator, site),
                .circle => |circle| if (!circle.cancelled) try self.circleEvent(allocator, circle),
            }
        }
    }

    pub fn addCircleEvent(self: *@This(), arc_idx: usize) !void {
        log.info("Adding circle event arc_idx = {d}, bounds of beachline are 0..{d}", .{ arc_idx, self.beachline.arcs.items.len });
        const beach_items = self.beachline.arcs.items;
        if (beach_items.len < 3) return;

        const arc = self.beachline.getArcPtr(arc_idx);

        const lo, const hi = get_focii_range: {
            if (arc.boundary_left == null) {
                log.info("Hit Start Cond", .{});
                break :get_focii_range .{ arc_idx, 3 };
            } else if (arc.boundary_right == null) {
                log.info("Hit End Cond", .{});
                break :get_focii_range .{ arc_idx - 2, beach_items.len };
            } else {
                log.info("Hit Middle Cond", .{});
                break :get_focii_range .{ arc_idx - 1, arc_idx + 2 };
            }
        };
        log.info("Range: {d}..{d}", .{ lo, hi });
        const focii = beach_items[lo..hi];

        if (circleFromPoints(focii[0].focus, focii[1].focus, focii[2].focus)) |circle| {
            self.beachline.getArcPtr(arc_idx).event = try self.event_queue.addTracking(.{
                .circle = .{
                    .center = circle.center,
                    .radius = circle.radius,
                    .arc = arc,
                },
            });
        } else |err| {
            log.err("Unable to generate circle {{{d}, {d}}}, {{{d}, {d}}}, {{{d}, {d}}}: {any}", .{
                focii[0].focus.x,
                focii[0].focus.y,
                focii[1].focus.x,
                focii[1].focus.y,
                focii[2].focus.x,
                focii[2].focus.y,
                err,
            });
        }
    }
    /// The only time an arc can appear on the beachline is during a site event
    /// Each event splits one arc in two
    pub fn siteEvent(self: *@This(), allocator: std.mem.Allocator, focus: std.meta.TagPayload(EventQueue.Event, .site)) !void {
        log.info("Beachline {any}", .{self.beachline});
        const beach_len = self.beachline.arcs.items.len;
        // Boundary Conditions
        // Locate the existing arc (if any) that is above the new site
        const idx, const yval =
            self.beachline.findOverlappingArcInfo(focus);

        if (self.beachline.getArcPtrOrNull(idx)) |old_arc| {
            // Invalidate the potential circle event in the event queue
            if (old_arc.event) |event| event.circle.cancelled = true;

            // Do some vector math
            const edge_start = rl.Vector2.init(focus.x, yval);
            const focus_offset: rl.Vector2 = .init(focus.x - old_arc.focus.x, focus.y - old_arc.focus.y);
            const edge_direction = rl.Vector2.init(focus_offset.y, -focus_offset.x).normalize();

            // Create initial edge data for the mesh
            const start_vertex = try self.mesh.addVertex(edge_start);
            const edge = try self.mesh.rawAddEdge(start_vertex, null, {}, true);

            // Dark magic tm.
            const should_not_have_left_boundary =
                (beach_len != 1 and old_arc.boundary_left == null) or
                (beach_len == 1 and focus.x < old_arc.focus.x);
            const should_not_have_right_boundary =
                (beach_len != 1 and old_arc.boundary_right == null) or
                (beach_len == 1 and focus.x > old_arc.focus.x);

            const inner_boundary_left = if (should_not_have_left_boundary)
                null
            else
                try self.beachline.createBoundary(.{
                    .mesh_hedge = edge.half,
                    .direction = edge_direction,
                    .start = edge_start,
                });

            const inner_boundary_right = if (should_not_have_right_boundary)
                null
            else
                try self.beachline.createBoundary(.{
                    .mesh_hedge = edge.half.twin,
                    .direction = edge_direction.negate(),
                    .start = edge_start,
                });

            const new_arcs: [3]Beachline.Arc = .{
                .{ // Left
                    .focus = old_arc.focus,
                    .event = null,
                    .boundary_left = old_arc.boundary_left,
                    .boundary_right = inner_boundary_left,
                },
                .{ // Middle
                    .focus = focus,
                    .event = null,
                    .boundary_left = inner_boundary_left,
                    .boundary_right = inner_boundary_right,
                },
                .{ // Right
                    .focus = old_arc.focus,
                    .event = null,
                    .boundary_left = inner_boundary_right,
                    .boundary_right = old_arc.boundary_right,
                },
            };

            const replacement_range: []const Beachline.Arc = if (should_not_have_left_boundary)
                new_arcs[1..3]
            else if (should_not_have_right_boundary)
                new_arcs[0..2]
            else
                &new_arcs;

            log.info("Replacing arc {{{d}, {d}}} at {d} with [{d}] -> {any}", .{
                focus.x,
                focus.y,
                idx,
                replacement_range.len,
                replacement_range,
            });

            try self.beachline.arcs.replaceRange(
                allocator,
                idx,
                1,
                replacement_range,
            );

            if (should_not_have_right_boundary) try self.addCircleEvent(
                idx -| 1,
            );
            if (should_not_have_left_boundary) try self.addCircleEvent(
                if (idx + 1 == beach_len) idx else idx + 1,
            );
        } else {
            @branchHint(.cold);
            try self.beachline.arcs.append(allocator, .{
                .focus = focus,
                .event = null,
                .boundary_left = null,
                .boundary_right = null,
            });
        }
        log.info("New Beachline {any}", .{self.beachline});
    }

    /// Circle events correspond to Voronoi vertices, false alarm events are never processed
    pub fn circleEvent(self: *@This(), allocator: std.mem.Allocator, inner: std.meta.TagPayload(EventQueue.Event, .circle)) !void {
        assert(inner.cancelled == false);
        assert(inner.cancelled == false);
        const arc: *Beachline.Arc = inner.arc;
        assert(arc.boundary_left != null and arc.boundary_right != null);
        const intersection_point = inner.center;
        _ = intersection_point; // autofix
        // inner: struct {
        //     center: rl.vector2,
        //     radius: f32,
        // },
        log.info("Running circle event", .{});
        _ = allocator; // autofix
        _ = self; // autofix
        // arc.boundary_left.?.mesh_hedge.();

        // Add the vertex to the corresponding edge in the mesh

        // Delete the arc and its circle events in the queue
        // Create a new edge in the mesh
        // Check the new triplets in the beachline for circle events
    }

    /// Computes the intersection of two parabolas given the directrix
    /// y - position of the directrix (sweepline)
    /// f1 - Focus of first parabola
    /// f2 - Focus of second parabola
    /// result - Intersection x-coordinate
    pub fn beachlineParabolaIntersection(directrix: f32, focusA: rl.Vector2, focusB: rl.Vector2) f32 {
        const fyDiff = focusA.y - focusB.y;
        if (fyDiff == 0) return (focusA.x + focusB.x) / 2;
        const fxDiff = focusA.x - focusB.x;
        const bAmd = focusA.y - directrix;
        const bBmd = focusB.y - directrix;
        const h1 = (-focusA.x * bBmd + focusB.x * bAmd) / fyDiff;
        const h2 = std.math.sqrt(bAmd * bBmd * (fxDiff * fxDiff + fyDiff * fyDiff)) / fyDiff;
        return h1 + h2;
    }
};

pub fn main() anyerror!void {
    const screenWidth = 800;
    const screenHeight = 600;

    var gpa: std.heap.GeneralPurposeAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    rl.initWindow(screenWidth, screenHeight, "Voronoi Cells");
    defer rl.closeWindow();

    rl.setTargetFPS(165);

    var points: std.ArrayListUnmanaged(rl.Vector2) = .empty;
    defer points.deinit(allocator);
    for (0..100) |_| {
        const padding = 100;
        const point = getRandomVector2(
            0 + padding / 2,
            screenWidth - padding / 2,
            0 + padding / 2,
            screenHeight - padding / 2,
        );

        try points.append(allocator, point);
    }

    // Voronoi
    var voronoi: Voronoi = .init(allocator);
    defer voronoi.deinit(allocator);

    try voronoi.generate(allocator, points.items, null);
    // while (voronoi.event_queue.count() > 0) try voronoi.processEvent(allocator);

    while (!rl.windowShouldClose()) {
        rl.beginDrawing();
        defer rl.endDrawing();

        rl.clearBackground(rl.Color.white);

        for (points.items) |point| {
            rl.drawCircleV(point, 4, rl.Color.gold);
        }

        rl.drawRectangleLinesEx(voronoi.clip_rect, 2, rl.Color.blue);
    }
}

fn drawParabola(focus: rl.Vector2, directrix: f32, minX: f32, maxX: f32, maxY: f32, color: rl.Color) void {
    _ = directrix; // autofix
    _ = minX; // autofix
    _ = maxX; // autofix
    _ = maxY; // autofix
    _ = color; // autofix
    const arc: Voronoi.Arc = .{
        .focus = focus,
        .squeeze_event = null,
    };
    _ = arc; // autofix

    const curvePts: [50]rl.Vector2 = undefined;
    _ = curvePts; // autofix
}

fn getRandomVector2(min_x: i32, max_x: i32, min_y: i32, max_y: i32) rl.Vector2 {
    return .init(
        @floatFromInt(rl.getRandomValue(min_x, max_x)),
        @floatFromInt(rl.getRandomValue(min_y, max_y)),
    );
}

const Complex = std.math.Complex(f32);
const Circle = struct { center: rl.Vector2, radius: f32 };
fn circleFromPoints(a: rl.Vector2, b: rl.Vector2, c: rl.Vector2) !Circle {
    const eq = std.meta.eql;
    const complex = std.math.complex;
    const z1: Complex = .init(a.x, a.y);
    const z2: Complex = .init(b.x, b.y);
    const z3: Complex = .init(c.x, c.y);
    if (eq(z1, z2) or eq(z2, z3) or eq(z1, z3))
        return error.DuplicatePoints;

    const z12 = z2.sub(z1);
    const w = z3.sub(z1).div(z12);

    // You should change 0 to a small tolerance for floating point comparisons
    if (std.math.approxEqAbs(f32, w.im, 0, std.math.floatEps(f32) * 5)) return error.PointsAreCollinear;

    const center: Complex = z12
        .mul(w.sub(.{ .re = complex.abs(w) * complex.abs(w), .im = 0 }))
        .div(w.sub(w.conjugate()))
        .add(z1); // Simplified denominator
    const radius: f32 = complex.abs(z1.sub(center));

    return .{ .center = .init(center.re, center.im), .radius = radius };
}
