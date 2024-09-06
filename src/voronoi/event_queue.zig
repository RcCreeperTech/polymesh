const std = @import("std");
const rl = @import("raylib");

const assert = std.debug.assert;

const Voronoi = @This();
const log = std.log.scoped(.voronoi_event_queue);

const EventQueue = @This();

pub const Event = union(enum) {
    const PriorityQueue = std.PriorityQueue(*@This(), void, lessThan);
    site: rl.Vector2,
    circle: struct {
        cancelled: bool = false,
        center: rl.Vector2,
        radius: f32,
        arc_idx: usize,
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

    const start_index = self.queue.items.len -| 1;
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
