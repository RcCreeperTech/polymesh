const std = @import("std");

const rl = @import("raylib");
const mesh = @import("mesh.zig");
const Voronoi = @import("voronoi.zig");

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
    for (0..10) |_| {
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
