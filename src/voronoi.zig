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

debug_last_event: ?*EventQueue.Event,

pub fn init(allocator: std.mem.Allocator) @This() {
    return .{
        .directrix = 0,
        .result_mesh = .empty,
        .beachline = .init(allocator),
        .event_queue = .init(allocator),
        .clip_rect = undefined, // FIXME:
        .debug_last_event = null,
    };
}

pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
    self.beachline.deinit(allocator);
    self.result_mesh.deinit(allocator);
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
            _ = try self.event_queue.add(
                .{ .site = point },
            );
        }
    }
    // TODO: Resize the clip rect

    while (self.event_queue.count() > 0) try self.processEvent(allocator);

    // TODO: Terminate remaining edges on the bounding box
}

pub fn debugDraw(self: *@This(), screen_w: f32) void {
    // Render The Beachline
    for (self.beachline.arc_list.items) |arc| {
        rl.drawCircleV(arc.parabola.focus, 5, rl.Color.blue);
        arc.parabola.draw(self.directrix, rl.Color.pink);
        arc.parabola.draw(self.directrix, rl.Color.pink);
    }
    // Render the last event
    if (self.debug_last_event) |ev|
        switch (ev.*) {
            .site => |site| {
                rl.drawCircleV(site, 7, rl.Color.red);
            },
            .circle => |circle| {
                rl.drawCircleLinesV(circle.center, circle.radius, rl.Color.green);
            },
        };
    // Render the sweepline
    rl.drawLineEx(
        rl.Vector2.init(0, self.directrix),
        rl.Vector2.init(screen_w, self.directrix),
        4,
        rl.Color.beige,
    );
}

pub fn processEvent(self: *@This(), allocator: std.mem.Allocator) !void {
    if (self.event_queue.removeOrNull()) |event| {
        log.info("Processing Event {any} ", .{event});
        self.debug_last_event = event;
        self.directrix = event.sortHandle();
        switch (event.*) {
            .site => |site| try self.siteEvent(allocator, site),
            .circle => |circle| if (!circle.cancelled) try self.circleEvent(circle),
        }
    }
}

pub fn addCircleEvent(self: *@This(), maybe_arc: ?*Beachline.Arc) !void {
    if (maybe_arc) |arc| {
        if (arc.arc_left == null or arc.arc_right == null) return;
        log.info("Trying to add circle event for arc {}, left = {} right = {}", .{ arc, arc.arc_left.?, arc.arc_right.? });

        if (circleFromPoints(arc.arc_left.?.parabola.focus, arc.parabola.focus, arc.arc_right.?.parabola.focus)) |circle| {
            log.info("Adding circle event, circle = {any}", .{circle});

            arc.event = try self.event_queue.add(.{
                .circle = .{
                    .center = circle.center,
                    .radius = circle.radius,
                    .arc = arc,
                },
            });
        } else |err| {
            log.err("Unable to generate circle {{{d}, {d}}}, {{{d}, {d}}}, {{{d}, {d}}}: {any}", .{
                arc.arc_left.?.parabola.focus.x,
                arc.arc_left.?.parabola.focus.y,
                arc.parabola.focus.x,
                arc.parabola.focus.y,
                arc.arc_right.?.parabola.focus.x,
                arc.arc_right.?.parabola.focus.y,
                err,
            });
        }
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
    if (self.beachline.findOverlappingArcInfoOrNull(focus)) |info| {
        const old_arc_idx, const old_arc_focus, const yval = info;
        // Do some vector math
        const edge_start = rl.Vector2.init(focus.x, yval);
        const focus_offset = focus.subtract(old_arc_focus);
        const edge_direction = rl.Vector2.init(focus_offset.y, -focus_offset.x).normalize();
        // Create initial edge data for the mesh
        const start_vertex = try self.result_mesh.addVertex(allocator, edge_start);
        // Insert the new Arc
        const new_arc = try self.beachline.insertArc(
            allocator,
            old_arc_idx,
            focus,
            .{
                .edge = try self.result_mesh.rawAddEdge(start_vertex, null, {}, true),
                .direction = edge_direction,
                .start = edge_start,
            },
            .{
                .edge = try self.result_mesh.rawAddEdge(start_vertex, null, {}, true),
                .direction = edge_direction.negate(),
                .start = edge_start,
            },
        );

        log.info("Beachline {any}", .{self.beachline});
        // Check for circle events to the left and right
        try self.addCircleEvent(new_arc.arc_left);
        try self.addCircleEvent(new_arc.arc_right);
    } else {
        _ = try self.beachline.insertArc(allocator, 0, focus, null, null);
        log.info("Beachline {any}", .{self.beachline});
    }
}

/// Circle events correspond to Voronoi vertices, false alarm events are never processed
pub fn circleEvent(self: *@This(), inner: std.meta.TagPayload(EventQueue.Event, .circle)) !void {
    assert(inner.cancelled == false);
    assert(inner.cancelled == false);
    const arc: *Beachline.Arc = inner.arc;
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
    const left = arc.arc_left.?;
    const right = arc.arc_right.?;
    // Remove the old arc and insert a boundary
    const offset = right.parabola.focus.subtract(left.parabola.focus);
    try self.beachline.removeArc(
        arc,
        .{
            .edge = try self.result_mesh.rawAddEdge(vertex, null, {}, true),
            .direction = rl.Vector2.init(offset.y, -offset.x).normalize(),
            .start = intersection_point,
        },
    );
    // Check the new triplets in the beachline for circle events
    try self.addCircleEvent(left);
    try self.addCircleEvent(right);
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
