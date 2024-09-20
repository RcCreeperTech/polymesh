const std = @import("std");
const rl = @import("raylib");
const poly = @import("polymesh");

const Mesh = poly.PolyMesh(VertData, void, [*:0]const u8);
const Delunator = poly.Delunator(Mesh);

const VertData = struct {
    pos: rl.Vector2,
    name: [*:0]const u8,

    pub fn format(
        self: @This(),
        comptime _: []const u8,
        _: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        return try writer.print("{{ .pos = .{{{d}, {d}}}, .name={s} }}", .{ self.pos.x, self.pos.y, self.name });
    }
};

pub fn main() !void {
    const screenWidth = 1600;
    const screenHeight = 900;

    const stdio = std.io.getStdOut().writer();

    var gpa: std.heap.GeneralPurposeAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    rl.initWindow(screenWidth, screenHeight, "Basic");
    defer rl.closeWindow();

    rl.setTargetFPS(60);
    rl.setWindowState(.{
        .msaa_4x_hint = true,
        .window_resizable = true,
    });

    var mesh: Mesh = .empty;
    defer mesh.deinit(allocator);

    const dd = try Delunator.init(&mesh);
    _ = dd; // autofix

    const v_names = [_][*:0]const u8{ "Va", "Vb", "Vc", "Vd", "Ve", "Vf" };
    const center: rl.Vector2 = .init(screenWidth / 2, screenHeight / 2);
    var shift: rl.Vector2 = .init(0, 200);

    var handles: std.ArrayListUnmanaged(Mesh.Vertex.Handle) = .empty;
    defer handles.deinit(allocator);
    try handles.append(allocator, try mesh.addVertex(allocator, .{ .pos = center, .name = "Vcenter" }));
    for (0..6) |i| {
        try handles.append(allocator, try mesh.addVertex(allocator, .{
            .pos = center.add(shift),
            .name = v_names[i],
        }));
        shift = shift.rotate(std.math.tau / 6.0);
    }

    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[1], {});
    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[2], {});
    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[3], {});
    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[4], {});
    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[5], {});
    // _ = try mesh.addEdge(allocator, handles.items[0], handles.items[6], {});

    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[1], handles.items[2] }, "A");
    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[2], handles.items[3] }, "B");
    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[3], handles.items[4] }, "C");
    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[4], handles.items[5] }, "D");
    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[5], handles.items[6] }, "E");
    _ = try mesh.addFace(allocator, &.{ handles.items[0], handles.items[6], handles.items[1] }, "F");

    try mesh.removeVertex(allocator, handles.items[0]);

    _ = try mesh.addFace(allocator, &.{
        handles.items[6],
        handles.items[5],
        handles.items[4],
        handles.items[3],
        handles.items[2],
        handles.items[1],
    }, "G");

    const hedges = &mesh.pools.half_edge;
    const edges = &mesh.pools.edge;
    const verts = &mesh.pools.vertex;
    const faces = &mesh.pools.face;

    try mesh.dumpDebugInfo(stdio, false);

    // try handles.items[2].free(allocator, vert_pool);

    while (!rl.windowShouldClose()) {
        rl.beginDrawing();
        defer rl.endDrawing();

        rl.clearBackground(rl.Color.fromInt(0x181818ff));

        {
            var it = edges.iterator();
            while (it.next()) |i| {
                drawEdge(&mesh, .{ .inner = i });
            }
        }

        {
            var it = hedges.iterator();
            while (it.next()) |i| {
                drawHalfEdge(&mesh, .{ .inner = i });
            }
        }

        {
            var it = verts.iterator();
            while (it.next()) |i| {
                const vert: Mesh.Vertex.Handle = .from(i);
                const vert_ptr = vert.deref(verts);
                const pos = vert_ptr.data.pos;
                const name = vert_ptr.data.name;
                rl.drawCircleV(pos, 5, rl.Color.white);
                rl.drawText(
                    name,
                    @intFromFloat(pos.x + 10),
                    @intFromFloat(pos.y + 10),
                    18,
                    .alpha(.yellow, 0.7),
                );
            }
        }

        {
            var it = faces.iterator();
            while (it.next()) |i| {
                const face: Mesh.Face.Handle = .from(i);
                var centroid: rl.Vector2 = .zero();

                const start = face.half(faces); // .twin(hedges);
                var current = start;
                var valence: f32 = 0;
                var size: f32 = 4;
                while (true) {
                    const Vp = current.origin(hedges).deref(verts).data.pos;
                    centroid = centroid.add(Vp);
                    rl.drawCircleV(Vp, size, rl.Color.alpha(.dark_green, 0.4));
                    size += 4;
                    valence += 1;

                    current = current.next(hedges);
                    if (current.inner == start.inner) break;
                }

                centroid = centroid.scale(1.0 / valence);
                // rl.drawCircleV(centroid, 5, rl.Color.dark_green);
                rl.drawText(
                    face.deref(faces).data,
                    @intFromFloat(centroid.x - 4),
                    @intFromFloat(centroid.y - 8),
                    18,
                    .dark_green,
                );
            }
        }
    }
}

const he_dist = 8;
fn drawEdge(mesh: *Mesh, edge: Mesh.Edge.Handle) void {
    const half = edge.half(&mesh.pools.edge);

    const origin = half
        .origin(&mesh.pools.half_edge)
        .deref(&mesh.pools.vertex);
    const destination = half
        .destination(&mesh.pools.half_edge)
        .deref(&mesh.pools.vertex);

    const start, const end = .{ origin.data.pos, destination.data.pos };
    const direction = end.subtract(start).normalize();

    // CCW rotation (to the left)
    const perp_shift = rl.Vector2
        .init(-direction.y, direction.x)
        .scale(he_dist);
    const midpoint = start.add(end).scale(0.5);

    // Render the link
    rl.drawLineEx(
        midpoint,
        midpoint.add(perp_shift),
        2,
        rl.Color.gold,
    );

    // Render the edge
    rl.drawLineEx(start, end, 2, rl.Color.ray_white);
}

fn drawHalfEdge(mesh: *Mesh, half: Mesh.HalfEdge.Handle) void {
    const origin = half
        .origin(&mesh.pools.half_edge)
        .deref(&mesh.pools.vertex);
    const destination = half
        .destination(&mesh.pools.half_edge)
        .deref(&mesh.pools.vertex);
    const start, const end = .{ origin.data.pos, destination.data.pos };

    const pad = 24;

    const direction = end.subtract(start).normalize();
    // CCW rotation (to the left)
    const perp_shift = rl.Vector2
        .init(-direction.y, direction.x)
        .scale(he_dist);
    const parallel_shift = direction
        .scale(pad);

    // Render the half edge
    rl.drawLineEx(
        start.add(perp_shift).add(parallel_shift),
        end.add(perp_shift).subtract(parallel_shift),
        2,
        rl.Color.red,
    );

    var tip = end.subtract(parallel_shift).add(perp_shift);
    // Render the half edge arrow
    rl.drawTriangle(
        tip,
        tip.add(perp_shift),
        tip.add(parallel_shift.scale(0.5)),
        rl.Color.red,
    );
}
