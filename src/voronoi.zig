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

pub fn init(allocator: std.mem.Allocator) @This() {
    return .{
        .directrix = 0,
        .result_mesh = .init(allocator),
        .beachline = .init(allocator),
        .event_queue = .init(allocator),
        .clip_rect = undefined, // FIXME:
    };
}

pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
    self.beachline.deinit(allocator);
    self.result_mesh.deinit();
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
    const beach_items = self.beachline.arcs.items;
    if (beach_items.len < 3) return;

    const arc = self.beachline.getArcPtr(arc_idx);
    if (arc.boundary_left == null or arc.boundary_right == null) return;

    const lo, const hi = .{ arc_idx - 1, arc_idx + 2 };
    log.info("Range: {d}..{d}", .{ lo, hi });
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
pub fn siteEvent(self: *@This(), allocator: std.mem.Allocator, focus: std.meta.TagPayload(EventQueue.Event, .site)) !void {
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
        const edge = try self.result_mesh.rawAddEdge(start_vertex, null, {}, true);

        const inner_boundary_left = try self.beachline.createBoundary(.{
            .mesh_hedge = edge.half,
            .direction = edge_direction,
            .start = edge_start,
        });

        const inner_boundary_right = try self.beachline.createBoundary(.{
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
        log.info("Beachline {any}", .{self.beachline});

        try self.addCircleEvent(idx -| 1);
        try self.addCircleEvent(idx + 2);
    } else {
        @branchHint(.cold);
        try self.beachline.arcs.append(allocator, .{
            .focus = focus,
            .event = null,
            .boundary_left = null,
            .boundary_right = null,
        });
    }
    log.info("Beachline {any}", .{self.beachline});
}

/// Circle events correspond to Voronoi vertices, false alarm events are never processed
pub fn circleEvent(self: *@This(), allocator: std.mem.Allocator, inner: std.meta.TagPayload(EventQueue.Event, .circle)) !void {
    assert(inner.cancelled == false);
    assert(inner.cancelled == false);
    const arc: *Beachline.Arc = self.beachline.getArcPtr(inner.arc_idx);
    assert(arc.boundary_left != null and arc.boundary_right != null);
    const intersection_point = inner.center;
    _ = intersection_point; // autofix
    // inner: struct {
    //     center: rl.vector2,
    //     radius: f32,
    // },
    log.info("Running circle event", .{});
    _ = allocator; // autofix
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
