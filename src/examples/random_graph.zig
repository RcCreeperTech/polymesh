const std = @import("std");
const rl = @import("raylib");
const PolyMesh = @import("polymesh").PolyMesh;

const Mesh = PolyMesh(struct { pos: rl.Vector2, vel: rl.Vector2 }, f32, void);

pub fn main() anyerror!void {
    const screenWidth = 1600;
    const screenHeight = 900;

    var gpa: std.heap.GeneralPurposeAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    rl.initWindow(screenWidth, screenHeight, "Basic");
    defer rl.closeWindow();

    rl.setTargetFPS(60);
    rl.setWindowState(.{ .msaa_4x_hint = true });

    var mesh: Mesh = .empty;
    defer mesh.deinit(allocator);
    const hedges = &mesh.pools.half_edge;
    const edges = &mesh.pools.edge;
    const verts = &mesh.pools.vertex;
    const faces = &mesh.pools.face;
    _ = faces; // autofix

    var handles: std.ArrayListUnmanaged(Mesh.Vertex.Handle) = .empty;
    defer handles.deinit(allocator);

    var points: std.ArrayListUnmanaged(rl.Vector2) = .empty;
    defer points.deinit(allocator);
    for (0..30) |_| {
        const padding = 100;
        const pos = getRandomVector2(
            0 + padding / 2,
            screenWidth - padding / 2,
            0 + padding / 2,
            screenHeight - padding / 2,
        );
        const vel = getRandomVector2(0, 5, 0, 5);
        try handles.append(
            allocator,
            try mesh.addVertex(allocator, .{ .pos = pos, .vel = vel }),
        );
    }

    for (handles.items) |from| {
        const rand: usize = @intCast(rl.getRandomValue(0, @intCast(handles.items.len - 1)));
        if (rl.getRandomValue(0, 100) > 50) {
            for (0..1) |_| {
                const to = handles.items[rand];
                _ = mesh.addEdge(
                    allocator,
                    from,
                    to,
                    from.deref(verts).data.pos.distance(to.deref(verts).data.pos),
                ) catch |err| switch (err) {
                    error.DuplicateVertices, error.EdgeAlreadyExists => {},
                    else => return err,
                };
            }
        }
    }

    while (!rl.windowShouldClose()) {
        // tick
        var vert_iter = verts.iterator();
        // Behold the worlds worst physics simulation
        while (vert_iter.next()) |vert| {
            var new_pos = vert.data.pos.add(vert.data.vel);

            var oob = false;
            if (new_pos.x > screenWidth or new_pos.x < 0) {
                vert.data.vel.x *= -1;
                oob = true;
            }
            if (new_pos.y > screenHeight or new_pos.y < 0) {
                vert.data.vel.y *= -1;
                oob = true;
            }

            if (vert.half) |half| {
                const constraint = half.parent(hedges).deref(edges).data;

                const other = half.destination(hedges).deref(verts);
                const distance = other.data.pos.distance(vert.data.pos);
                const err = distance - constraint;
                if (err != 0) {
                    const correction = other.data.pos
                        .subtract(vert.data.pos)
                        .normalize()
                        .scale(err * 0.5);
                    new_pos = new_pos.add(correction);

                    const new_other_pos = other.data.pos.subtract(correction);
                    var ooob = false;
                    if (new_other_pos.x > screenWidth or new_other_pos.x < 0) {
                        other.data.vel.x *= -1;
                        ooob = true;
                    }
                    if (new_other_pos.y > screenHeight or new_other_pos.y < 0) {
                        other.data.vel.y *= -1;
                        ooob = true;
                    }
                    if (!ooob) other.data.pos = new_other_pos;
                }
            }

            if (!oob) vert.data.pos = new_pos;
        }

        rl.beginDrawing();
        defer rl.endDrawing();

        rl.clearBackground(rl.Color.fromInt(0x181818ff));

        var edge_iter = edges.iterator();
        while (edge_iter.next()) |edge| {
            const half = edge.half;
            const from = half.origin(hedges).deref(verts);
            const to = half.twin(hedges).origin(hedges).deref(verts);
            rl.drawLineEx(from.data.pos, to.data.pos, 2, rl.Color.ray_white);
        }

        vert_iter.reset();
        while (vert_iter.next()) |vert| {
            rl.drawCircleV(vert.data.pos, 5, rl.Color.red);
        }
    }
}

fn getRandomVector2(min_x: i32, max_x: i32, min_y: i32, max_y: i32) rl.Vector2 {
    return .init(
        @floatFromInt(rl.getRandomValue(min_x, max_x)),
        @floatFromInt(rl.getRandomValue(min_y, max_y)),
    );
}
