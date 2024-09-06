const std = @import("std");
const rl = @import("raylib");

const mesh = @import("../mesh.zig");
const EventQueue = @import("event_queue.zig");

const assert = std.debug.assert;
pub const PolyMesh = mesh.GenericPolyMesh(rl.Vector2, void, void);

const log = std.log.scoped(.voronoi_beachline);

pub const Arc = struct {
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
    if (self.arcs.items.len == 0) return .{ 0, 0 };

    const higher_idx = std.sort.lowerBound(Arc, self.arcs.items, point.x, Arc.compare);
    const lower_idx = higher_idx -| 1;

    const lower_arc = self.getArc(lower_idx);
    const ylow = lower_arc.evaluate(point.x, point.y);

    if (higher_idx != 0 and higher_idx != self.arcs.items.len) {
        const yhigh = self.getArc(higher_idx).evaluate(point.x, point.y);
        log.info("Resolving tie lower_idx = {d}, low = {d}, high = {d}", .{ lower_idx, ylow, yhigh });
        if (yhigh > ylow) {
            return .{ higher_idx, yhigh };
        } else {
            return .{ lower_idx, ylow };
        }
    }

    return .{ lower_idx, ylow };
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
