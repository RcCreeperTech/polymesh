const std = @import("std");
const rl = @import("raylib");

const mesh = @import("../polymesh.zig");
const EventQueue = @import("event_queue.zig");

const assert = std.debug.assert;
/// A null value corresponds to the point at infinity
/// Having multiple vertices with this value is undefined behavior
pub const PolyMesh = mesh.PolyMesh(?rl.Vector2, void, void);

const log = std.log.scoped(.voronoi_beachline);

pub const Parabola = struct {
    focus: rl.Vector2,

    pub fn evaluate(self: *const @This(), x: f32, directrix_y: f32) f32 {
        // a = focus.x
        // b = focus.y
        // c = directrix
        // distance to focus = distance to directrix
        // sqrt[(x - a)^2 + (y - b)^2] = |y - c|
        // (x - a)^2 + (y - b)^2 = (y - c)^2
        // (x - a)^2 + y^2 - 2yb +b^2 = y^2 - 2yc + c^2
        // (x - a)^2 + b^2 - c^2 = y^2 - 2yc - y^2 + 2yb
        // (x - a)^2 + b^2 - c^2 = - 2yc + 2yb
        // (x - a)^2 + b^2 - c^2 = (-2c + 2b)y
        // [(x - a)^2 + b^2 - c^2]/(-2c + 2b) = y
        // y = [(x - a)^2 + b^2 - c^2]/(-2c + 2b)
        const focus = self.focus;
        const diff_x = x - focus.x;
        var y = diff_x * diff_x + focus.y * focus.y - directrix_y * directrix_y;
        y /= -2 * directrix_y + 2 * focus.y;
        return y;
    }

    pub fn draw(self: *const @This(), directrix: f32, color: rl.Color) void {
        const bound = 400;
        const step: f32 = 15;
        if (std.math.approxEqAbs(f32, directrix, self.focus.y, std.math.floatEps(f32))) {
            rl.drawLineEx(self.focus, .{ .x = self.focus.x, .y = 0 }, 3, rl.Color.dark_purple);
            return;
        }
        self.drawEx(directrix, color, bound, step);
    }

    pub fn drawEx(self: *const @This(), directrix: f32, color: rl.Color, bound: f32, step_size: f32) void {
        const focus = self.focus;
        const p: Parabola = .{ .focus = focus };
        var x: f32 = p.focus.x - bound;
        var y: f32 = p.evaluate(x, directrix);
        while (x < focus.x + bound) : (x += step_size) {
            const next_y = p.evaluate(x + step_size, directrix);
            rl.drawLineEx(.{ .x = x, .y = y }, .{ .x = x + step_size, .y = next_y }, 4, color);
            y = next_y;
        }
    }

    pub fn format(
        self: @This(),
        comptime _: []const u8,
        _: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        return writer.print("{{{d}, {d}}}", .{ self.focus.x, self.focus.y });
    }
};

pub const Arc = struct {
    parabola: Parabola,
    event: ?*EventQueue.Event,
    boundary_left: ?*Boundary,
    boundary_right: ?*Boundary,
    arc_left: ?*Arc,
    arc_right: ?*Arc,

    pub const default: @This() = .{
        .focus = undefined,
        .event = null,
        .boundary_left = null,
        .boundary_right = null,
    };

    fn compare(x: f32, item: *@This()) std.math.Order {
        return std.math.order(item.parabola.focus.x, x);
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

        try writer.print("{}", .{self.parabola});
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
    edge: *PolyMesh.Edge,
};

/// maintains the arc ordering by x-coordinate
arc_list: std.ArrayListUnmanaged(*Arc),
arc_pool: std.heap.MemoryPool(Arc),
boundary_pool: std.heap.MemoryPool(Boundary),

pub fn init(allocator: std.mem.Allocator) @This() {
    return .{
        .arc_list = .empty,
        .arc_pool = .init(allocator),
        .boundary_pool = .init(allocator),
    };
}

pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
    self.arc_list.deinit(allocator);
    self.arc_pool.deinit();
    self.boundary_pool.deinit();
}

pub fn findOverlappingArcInfoOrNull(self: *@This(), point: rl.Vector2) ?struct {
    usize, // arc_idx
    rl.Vector2, // arc_focus
    f32, // yval
} {
    log.info("Searcing Point {{{d}, {d}}}", .{ point.x, point.y });
    if (self.arc_list.items.len == 0) return null;

    const higher_idx = std.sort.lowerBound(*Arc, self.arc_list.items, point.x, Arc.compare);
    const lower_idx = higher_idx -| 1;

    const lower_arc = self.getArcPtr(lower_idx);
    const ylow = lower_arc.parabola.evaluate(point.x, point.y);

    if (higher_idx != 0 and higher_idx != self.arc_list.items.len) {
        const higher_arc = self.getArcPtr(higher_idx);
        const yhigh = higher_arc.parabola.evaluate(point.x, point.y);
        log.info("Resolving tie lower_idx = {d}, low = {d}, high = {d}", .{ lower_idx, ylow, yhigh });
        if (yhigh > ylow) {
            return .{ higher_idx, higher_arc.parabola.focus, yhigh };
        } else {
            return .{ lower_idx, lower_arc.parabola.focus, ylow };
        }
    }

    return .{ lower_idx, lower_arc.parabola.focus, ylow };
}

pub fn getArcPtr(self: *const @This(), index: usize) *Arc {
    return self.arc_list.items[index];
}

pub fn getArcPtrOrNull(self: *const @This(), index: i32) ?*Arc {
    return if (index >= self.arc_list.items.len or index < 0)
        null
    else
        self.arc_list.items[@intCast(index)];
}

pub fn insertArc(
    self: *@This(),
    allocator: std.mem.Allocator,
    replacement_idx: usize,
    focus: rl.Vector2,
    boundary_left: ?Boundary,
    boundary_right: ?Boundary,
) !*Arc {
    const new_boundary_left = if (boundary_left) |boundary| try self.createBoundary(boundary) else null;
    const new_boundary_right = if (boundary_right) |boundary| try self.createBoundary(boundary) else null;
    if (self.arc_list.items.len != 0) {
        const arc = self.getArcPtr(replacement_idx);
        // Invalidate the potential circle event in the event queue
        if (arc.event) |event| event.circle.cancelled = true;

        const left_arc = self.getArcPtrOrNull(@as(i32, @intCast(replacement_idx)) - 1);
        const right_arc = self.getArcPtrOrNull(@as(i32, @intCast(replacement_idx + 1)));

        const new_arc_left = try self.arc_pool.create();
        const new_arc_middle = try self.arc_pool.create();
        const new_arc_right = try self.arc_pool.create();

        if (left_arc) |left| left.arc_right = new_arc_left;
        new_arc_left.* = .{
            .parabola = arc.parabola,
            .event = null,
            .boundary_left = arc.boundary_left,
            .boundary_right = new_boundary_left,
            .arc_left = left_arc,
            .arc_right = new_arc_middle,
        };

        new_arc_middle.* = .{
            .parabola = .{ .focus = focus },
            .event = null,
            .boundary_left = new_boundary_left,
            .boundary_right = new_boundary_right,
            .arc_left = new_arc_left,
            .arc_right = new_arc_right,
        };

        new_arc_right.* = .{
            .parabola = arc.parabola,
            .event = null,
            .boundary_left = new_boundary_right,
            .boundary_right = arc.boundary_right,
            .arc_left = new_arc_middle,
            .arc_right = right_arc,
        };
        if (right_arc) |right| right.arc_left = new_arc_right;

        const new_arcs = &.{ new_arc_left, new_arc_middle, new_arc_right };
        log.info("Inserting Arcs {} at index {d}", .{ new_arcs, replacement_idx });
        try self.arc_list.replaceRange(allocator, replacement_idx, 1, new_arcs);

        self.arc_pool.destroy(arc);
        return new_arc_middle;
    } else {
        @branchHint(.cold);
        const new_arc = try self.arc_pool.create();
        new_arc.* = .{
            .parabola = .{ .focus = focus },
            .event = null,
            .boundary_left = new_boundary_left,
            .boundary_right = new_boundary_right,
            .arc_left = null,
            .arc_right = null,
        };
        try self.arc_list.append(allocator, new_arc);
        return new_arc;
    }
}

pub fn removeArc(self: *@This(), arc: *Arc, merge_boundary: Boundary) !void {
    // PERF: Can I get rid of this or speed it up?
    const removal_idx = found_idx: for (self.arc_list.items, 0..) |a, i| {
        if (a == arc) break :found_idx i;
    } else unreachable;

    if (arc.arc_left != null and arc.arc_right != null) {
        const left = arc.arc_left.?;
        const right = arc.arc_right.?;
        // Remove potential events
        if (left.event) |event| event.circle.cancelled = true;
        if (right.event) |event| event.circle.cancelled = true;

        const new_boundary = try self.createBoundary(merge_boundary);
        // Stitch up the gap
        left.boundary_right = new_boundary;
        right.boundary_left = new_boundary;

        left.arc_right = right;
        right.arc_left = left;
        // Destroy old boundaries
        self.boundary_pool.destroy(arc.boundary_left.?);
        self.boundary_pool.destroy(arc.boundary_right.?);
    }
    // Destroy the arc
    _ = self.arc_list.orderedRemove(removal_idx);
    self.arc_pool.destroy(arc);
}

fn createBoundary(self: *@This(), boundary: Boundary) !*Boundary {
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
    for (self.arc_list.items, 0..) |arc, i| {
        try writer.print("{any}", .{arc});
        if (i != self.arc_list.items.len - 1) try writer.print(",", .{});
    }
    return writer.print("]", .{});
}
