const std = @import("std");

const rl = @import("raylib");
const mesh = @import("mesh.zig");
const Voronoi = @import("voronoi.zig");
const Parabola = @import("voronoi/beachline.zig").Parabola;

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
    const screenWidth = 1600;
    const screenHeight = 900;

    var gpa: std.heap.GeneralPurposeAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    rl.initWindow(screenWidth, screenHeight, "Voronoi Cells");
    defer rl.closeWindow();

    rl.setTargetFPS(165);

    var points: std.ArrayListUnmanaged(rl.Vector2) = .empty;
    defer points.deinit(allocator);
    for (0..30) |_| {
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

    voronoi.point_at_infinity = try voronoi.result_mesh.addVertex(null);
    voronoi.clip_rect = gen_clip_rect: {
        const aabb_padding = 20;
        const infinity = std.math.inf(f32);
        var min, var max = .{
            rl.Vector2.init(infinity, infinity),
            rl.Vector2.init(-infinity, -infinity),
        };

        for (points.items) |point| {
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
    };

    for (points.items) |point| {
        if (rl.checkCollisionPointRec(point, voronoi.clip_rect)) {
            _ = try voronoi.event_queue.add(
                .{ .site = point },
            );
        }
    }
    // TODO: Resize the clip rect

    while (voronoi.event_queue.count() > 0) try voronoi.processEvent(allocator);

    // TODO: Terminate remaining edges on the bounding box

    // while (voronoi.event_queue.count() > 0) try voronoi.processEvent(allocator);

    while (!rl.windowShouldClose()) {
        rl.beginDrawing();
        defer rl.endDrawing();

        const mouse_y: f32 = @floatFromInt(rl.getMouseY());

        rl.clearBackground(rl.Color.white);

        for (points.items) |point| {
            rl.drawCircleV(point, 4, rl.Color.red);
            const p: Parabola = .{ .focus = point };
            p.draw(mouse_y, rl.Color.pink);
        }

        rl.drawLineEx(
            rl.Vector2.init(0, mouse_y),
            rl.Vector2.init(screenWidth, mouse_y),
            4,
            rl.Color.gold,
        );

        var it = voronoi.result_mesh.vertex_list.first;
        while (it) |node| : (it = node.next) {
            const vert = node.data;
            if (vert.data) |point| {
                rl.drawCircleV(point, 4, rl.Color.green);
                if (vert.half.?.twin.origin.data) |b| {
                    rl.drawLineV(point, b, rl.Color.blue);
                }
            }
        }
    }
}

fn getRandomVector2(min_x: i32, max_x: i32, min_y: i32, max_y: i32) rl.Vector2 {
    return .init(
        @floatFromInt(rl.getRandomValue(min_x, max_x)),
        @floatFromInt(rl.getRandomValue(min_y, max_y)),
    );
}
