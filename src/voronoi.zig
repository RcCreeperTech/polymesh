const std = @import("std");
const rl = @import("raylib");

const EventQueue = @import("voronoi/event_queue.zig");
const Beachline = @import("voronoi/beachline.zig");

const assert = std.debug.assert;

const log = std.log.scoped(.voronoi);

const Voronoi = @This();

directrix: f32,
result_mesh: Beachline.PolyMesh,
beachline: Beachline,
event_queue: EventQueue,
clip_rect: rl.Rectangle,
point_at_infinity: *Beachline.PolyMesh.Vertex,

pub fn init(allocator: std.mem.Allocator) @This() {
    return .{
        .directrix = 0,
        .result_mesh = .init(allocator),
        .beachline = .init(allocator),
        .event_queue = .init(allocator),
        .clip_rect = undefined, // FIXME:
        .point_at_infinity = undefined,
    };
}

pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
    self.beachline.deinit(allocator);
    self.result_mesh.deinit();
    self.event_queue.deinit();
}

pub fn generate(self: *@This(), allocator: std.mem.Allocator, points: []rl.Vector2, optional_clip_rect: ?rl.Rectangle) !void {
    self.point_at_infinity = try self.result_mesh.addVertex(null);
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
            .circle => |circle| if (!circle.cancelled) try self.circleEvent(circle),
        }
    }
}

pub fn addCircleEvent(self: *@This(), arc_idx: usize) !void {
    const beach_items = self.beachline.arcs.items;
    if (beach_items.len < 3) return;
    if (arc_idx == 0 or arc_idx == beach_items.len - 1) return;

    const lo, const hi = .{ arc_idx - 1, arc_idx + 2 };
    log.info("Attempting to add Circle event with range = {d}..{d}", .{ lo, hi });
    const focii = beach_items[lo..hi];

    if (circleFromPoints(focii[0].focus, focii[1].focus, focii[2].focus)) |circle| {
        log.info("Adding circle event arc_idx = {d}, beach_len = {d}, circle = {any}", .{ arc_idx, self.beachline.arcs.items.len, circle });
        self.beachline.getArcPtr(arc_idx).event = try self.event_queue.addTracking(.{
            .circle = .{
                .center = circle.center,
                .radius = circle.radius,
                .arc_idx = arc_idx,
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
pub fn siteEvent(
    self: *@This(),
    allocator: std.mem.Allocator,
    focus: std.meta.TagPayload(EventQueue.Event, .site),
) !void {
    log.info("Beachline {any}", .{self.beachline});
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
        const start_vertex = try self.result_mesh.addVertex(edge_start);

        const inner_boundary_left = try self.beachline.createBoundary(.{
            .edge = try self.result_mesh.rawAddEdge(start_vertex, null, {}, true),
            .direction = edge_direction,
            .start = edge_start,
        });

        const inner_boundary_right = try self.beachline.createBoundary(.{
            .edge = try self.result_mesh.rawAddEdge(start_vertex, null, {}, true),
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

        const replacement_range: []const Beachline.Arc = &new_arcs;

        log.info("Replacing arc at index {d} with {any}", .{
            idx,
            replacement_range,
        });

        try self.beachline.arcs.replaceRange(
            allocator,
            idx,
            1,
            replacement_range,
        );
        // PERF: Currently we need to fixup all the circle events that refer to all
        //       arcs to the right of the added ones because they are holding indices
        //       into a flat list of arcs :(
        for (self.event_queue.queue.items) |item| switch (item.*) {
            .circle => |circle| {
                if (circle.arc_idx > idx) item.circle.arc_idx += 2;
            },
            else => {},
        };
        log.info("Beachline {any}", .{self.beachline});

        const new_arc_idx = idx + 1; // beacause we always insert "old new old" on top of "old"
        try self.addCircleEvent(new_arc_idx - 1);
        try self.addCircleEvent(new_arc_idx + 1);
    } else {
        @branchHint(.cold);
        try self.beachline.arcs.append(allocator, .{
            .focus = focus,
            .event = null,
            .boundary_left = try self.beachline.createBoundary(.{
                .edge = try self.result_mesh.rawAddEdge(self.point_at_infinity, null, {}, true),
                .direction = .init(0, 1),
                .start = .zero(),
            }),
            .boundary_right = try self.beachline.createBoundary(.{
                .edge = try self.result_mesh.rawAddEdge(self.point_at_infinity, null, {}, true),
                .direction = .init(0, 1),
                .start = .zero(),
            }),
        });
    }
    log.info("Beachline {any}", .{self.beachline});
}

/// Circle events correspond to Voronoi vertices, false alarm events are never processed
pub fn circleEvent(self: *@This(), inner: std.meta.TagPayload(EventQueue.Event, .circle)) !void {
    assert(inner.cancelled == false);
    assert(inner.cancelled == false);
    const arc: *Beachline.Arc = self.beachline.getArcPtr(inner.arc_idx);
    assert(arc.boundary_left != null and arc.boundary_right != null);
    const intersection_point = inner.center;
    log.info("Running circle event", .{});

    // Add the vertex to the corresponding edge in the mesh
    const vertex = try self.result_mesh.addVertex(intersection_point);
    // Link the vertex to the right edge
    vertex.half = arc.boundary_right.?.edge.half.twin;
    // Link the edges together
    arc.boundary_left.?.edge.half.twin.origin = vertex;
    arc.boundary_right.?.edge.half.twin.origin = vertex;
    // Get the left and right arcs
    const left = self.beachline.getArcPtr(inner.arc_idx - 1);
    const right = self.beachline.getArcPtr(inner.arc_idx + 1);
    // Create the new boundary to join the two arcs
    const offset = right.focus.subtract(left.focus);
    const new_boundary = try self.beachline.createBoundary(.{
        .edge = try self.result_mesh.rawAddEdge(vertex, null, {}, true),
        .direction = rl.Vector2.init(offset.y, -offset.x).normalize(),
        .start = intersection_point,
    });
    // Stitch up the gap before we invalidate our pointers
    left.boundary_right = new_boundary;
    right.boundary_left = new_boundary;
    // Destroy merged boundaries
    self.beachline.boundary_pool.destroy(arc.boundary_left.?);
    self.beachline.boundary_pool.destroy(arc.boundary_right.?);
    // Remove the arc and the duplicate arc to the right
    // NOTE: This operation invalidates arcPtr's and indexs above the removed arc
    _ = self.beachline.arcs.orderedRemove(inner.arc_idx);
    // PERF: Currently we need to fixup all the circle events that refer to all
    //       arcs to the right of the removed one because they are holding indices
    //       into a flat list of arcs
    for (self.event_queue.queue.items) |item| switch (item.*) {
        .circle => |circle| {
            if (circle.arc_idx > inner.arc_idx) item.circle.arc_idx += 1;
        },
        else => {},
    };
    // Check the new triplets in the beachline for circle events
    const left_idx = inner.arc_idx - 1; // NOTE: This should never underflow
    const right_idx = inner.arc_idx; // We did a remove on the old arc so this is now the correct index
    try self.addCircleEvent(left_idx);
    try self.addCircleEvent(right_idx);
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
