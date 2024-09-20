const std = @import("std");
const mem = std.mem;
const Allocator = mem.Allocator;

const SlotPool = @import("slot_pool.zig").SlotPool;

pub fn PolyMesh(comptime VertexData: type, comptime EdgeData: type, comptime FaceData: type) type {
    return struct {
        const Self = @This();

        pub const Vertex = struct {
            pub const Pool = SlotPool(@This(), .{});
            pub const Handle = struct {
                inner: Pool.Handle,
                pub fn from(h: Pool.Handle) @This() {
                    return .{ .inner = h };
                }

                pub fn alloc(allocator: Allocator, pool: *Vertex.Pool) !@This() {
                    return .{ .inner = try pool.alloc(allocator) };
                }
                pub fn free(self: *const @This(), allocator: Allocator, pool: *Vertex.Pool) !void {
                    return try pool.free(allocator, self.inner);
                }
                /// Convinience function for half edge traversal
                pub fn deref(self: *const @This(), pool: *Vertex.Pool) *Vertex {
                    return self.inner.deref(pool);
                }
                /// Convinience function for vertex traversal
                pub fn halfOrNull(self: *const @This(), pool: *Vertex.Pool) ?HalfEdge.Handle {
                    return self.inner.deref(pool).half;
                }

                pub fn format(self: @This(), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
                    return std.fmt.formatIntValue(@intFromEnum(self.inner), fmt, options, writer);
                }
            };
            half: ?HalfEdge.Handle, // The edge linked must always be outbound
            data: VertexData,
        };

        pub const Edge = struct {
            pub const Pool = SlotPool(@This(), .{});
            pub const Handle = struct {
                inner: Pool.Handle,
                pub fn from(h: Pool.Handle) @This() {
                    return .{ .inner = h };
                }
                pub fn alloc(allocator: Allocator, pool: *Edge.Pool) !@This() {
                    return .{ .inner = try pool.alloc(allocator) };
                }
                pub fn free(self: *const @This(), allocator: Allocator, pool: *Edge.Pool) !void {
                    return try pool.free(allocator, self.inner);
                }
                /// Convinience function for half edge traversal
                pub fn deref(self: *const @This(), pool: *Edge.Pool) *Edge {
                    return self.inner.deref(pool);
                }
                /// Convinience function for vertex traversal
                pub fn half(self: *const @This(), pool: *Edge.Pool) HalfEdge.Handle {
                    return self.inner.deref(pool).half;
                }

                pub fn format(self: @This(), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
                    return std.fmt.formatIntValue(@intFromEnum(self.inner), fmt, options, writer);
                }
            };
            half: HalfEdge.Handle,
            data: EdgeData,
        };

        pub const Face = struct {
            pub const Pool = SlotPool(@This(), .{});
            pub const Handle = struct {
                inner: Pool.Handle,
                pub fn from(h: Pool.Handle) @This() {
                    return .{ .inner = h };
                }
                pub fn alloc(allocator: Allocator, pool: *Face.Pool) !@This() {
                    return .{ .inner = try pool.alloc(allocator) };
                }
                pub fn free(self: *const @This(), allocator: Allocator, pool: *Face.Pool) !void {
                    return try pool.free(allocator, self.inner);
                }
                /// Convinience function for half edge traversal
                pub fn deref(self: *const @This(), pool: *Face.Pool) *Face {
                    return self.inner.deref(pool);
                }
                /// Convinience function for vertex traversal
                pub fn half(self: *const @This(), pool: *Face.Pool) HalfEdge.Handle {
                    return self.inner.deref(pool).half;
                }

                pub const VertexIterator = struct {
                    current: ?Vertex.Handle,
                    start: Vertex.Handle,
                    mesh: *Self,

                    pub fn next(self: *@This()) ?*Vertex {
                        std.debug.print("C: {?}\n", .{self.current});
                        if (self.current) |current| { // i > 0
                            const next_vert = current
                                .halfOrNull(&self.mesh.pools.vertex).?
                                .destination(&self.mesh.pools.half_edge);
                            if (next_vert.inner == current.inner) return null;

                            defer self.current = next_vert;
                            return current.deref(&self.mesh.pools.vertex);
                        } else { // i = 0
                            self.current = self.start;
                            return self.start.deref(&self.mesh.pools.vertex);
                        }
                    }

                    /// Reset the iterator to the initial index
                    pub fn reset(self: *@This()) void {
                        self.current_handle = self.start_handle;
                    }
                };

                pub fn vertexIter(self: *const @This(), mesh: *Self) VertexIterator {
                    return .{
                        .current = null,
                        .start = self
                            .half(&mesh.pools.face)
                            .origin(&mesh.pools.half_edge),
                        .mesh = mesh,
                    };
                }

                pub fn format(self: @This(), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
                    return std.fmt.formatIntValue(@intFromEnum(self.inner), fmt, options, writer);
                }
            };
            half: HalfEdge.Handle,
            data: FaceData,
        };

        pub const HalfEdge = struct {
            pub const Pool = SlotPool(@This(), .{});
            pub const Handle = struct {
                inner: Pool.Handle,
                pub fn from(h: Pool.Handle) @This() {
                    return .{ .inner = h };
                }
                pub fn alloc(allocator: Allocator, pool: *HalfEdge.Pool) !@This() {
                    return .{ .inner = try pool.alloc(allocator) };
                }
                pub fn free(self: *const @This(), allocator: Allocator, pool: *HalfEdge.Pool) !void {
                    return try pool.free(allocator, self.inner);
                }
                /// Convinience function for half edge traversal
                pub fn deref(self: *const @This(), pool: *HalfEdge.Pool) *HalfEdge {
                    return self.inner.deref(pool);
                }
                /// Convinience function for half edge traversal
                pub fn origin(self: *const @This(), pool: *HalfEdge.Pool) Vertex.Handle {
                    return self.inner.deref(pool).origin;
                }
                /// Convinience function for half edge traversal
                pub fn destination(self: *const @This(), pool: *HalfEdge.Pool) Vertex.Handle {
                    return self.twin(pool).origin(pool);
                }
                /// Convinience function for half edge traversal
                pub fn leftFaceOrNull(self: *const @This(), pool: *HalfEdge.Pool) ?Face.Handle {
                    return self.inner.deref(pool).left_face;
                }
                /// Convinience function for half edge traversal
                pub fn parent(self: *const @This(), pool: *HalfEdge.Pool) Edge.Handle {
                    return self.inner.deref(pool).parent;
                }
                /// Convinience function for half edge traversal
                pub fn twin(self: *const @This(), pool: *HalfEdge.Pool) HalfEdge.Handle {
                    return self.inner.deref(pool).twin;
                }
                /// Convinience function for half edge traversal
                pub fn next(self: *const @This(), pool: *HalfEdge.Pool) HalfEdge.Handle {
                    return self.inner.deref(pool).next;
                }
                /// Convinience function for half edge traversal
                pub fn prev(self: *const @This(), pool: *HalfEdge.Pool) HalfEdge.Handle {
                    return self.inner.deref(pool).prev;
                }

                pub fn format(self: @This(), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
                    return std.fmt.formatIntValue(@intFromEnum(self.inner), fmt, options, writer);
                }
            };
            origin: Vertex.Handle,
            left_face: ?Face.Handle,
            parent: Edge.Handle,

            twin: HalfEdge.Handle,
            next: HalfEdge.Handle,
            prev: HalfEdge.Handle,
        };

        verts: Vertex.Pool,
        edges: Edge.Pool,
        faces: Face.Pool,
        hedges: HalfEdge.Pool,

        pub const empty: Self = .{
            .verts = .empty,
            .faces = .empty,
            .edges = .empty,
            .hedges = .empty,
        };

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.verts.deinit(allocator);
            self.edges.deinit(allocator);
            self.faces.deinit(allocator);
            self.hedges.deinit(allocator);
        }

        pub fn dumpDebugInfo(self: *Self, writer: anytype, comptime show_dead: bool) !void {
            const hedges = &self.hedges;
            const edges = &self.edges;
            const verts = &self.verts;
            const faces = &self.faces;

            const table_names = .{ "Vertex", "Face", "Edge", "Half Edge" };
            const tables = .{ verts, faces, edges, hedges };

            inline for (tables, 0..) |table, tdx| {
                const easy_table = tdx != 3;

                const table_prelude = if (show_dead and easy_table)
                    \\
                    \\ Idx | Ded | Half | Data
                    \\-----|-----|------|------------------------
                else if (show_dead and !easy_table)
                    \\
                    \\ Idx | Ded | Twin | Next | Prev | Parent | Left Face | Origin
                    \\-----|-----|------|------|------|--------|-----------|--------
                else if (!show_dead and easy_table)
                    \\
                    \\ Idx | Half | Data
                    \\-----|------|------------------------
                else
                    \\
                    \\ Idx | Twin | Next | Prev | Parent | Left Face | Origin
                    \\-----|------|------|------|--------|-----------|--------
                    ;

                try writer.print(table_names[tdx] ++ " Table\n" ++ table_prelude ++ "\n", .{});

                if (show_dead) {
                    for (table.slot_list.items, 0..) |item, idx| {
                        if (easy_table) {
                            try writer.print(" {d: ^3} | {s: ^3} | {?: ^4} | {any}\n", .{
                                idx,
                                if (table.free_bitset.isSet(idx)) "X" else "",
                                item.half,
                                item.data,
                            });
                        } else {
                            try writer.print(" {d: ^3} | {s: ^3} | {d: ^4} | {d: ^4} | {d: ^4} | {d: ^6} | {?: ^9} | {d: ^6}\n", .{
                                idx,
                                if (table.free_bitset.isSet(idx)) "X" else "",
                                item.twin,
                                item.next,
                                item.prev,
                                item.parent,
                                item.left_face,
                                item.origin,
                            });
                        }
                    }
                } else {
                    for (table.slot_list.items, 0..) |item, idx| {
                        if (!table.free_bitset.isSet(idx)) {
                            if (easy_table) {
                                try writer.print(" {d: ^3} | {?: ^4} | {any}\n", .{
                                    idx,
                                    item.half,
                                    item.data,
                                });
                            } else {
                                try writer.print(" {d: ^3} | {d: ^4} | {d: ^4} | {d: ^4} | {d: ^6} | {?: ^9} | {d: ^6}\n", .{
                                    idx,
                                    item.twin,
                                    item.next,
                                    item.prev,
                                    item.parent,
                                    item.left_face,
                                    item.origin,
                                });
                            }
                        }
                    }
                }
                try writer.print("\n", .{});
            }
        }

        /// Given a vertex finds the first outbound half edge without a left face.
        /// Returns error in the case that there is not a half edge satisfying this property
        pub fn findFirstFreeOutboundHalfEdgeOrNull(self: *Self, vert: Vertex.Handle) !?HalfEdge.Handle {
            const hedges = &self.pools.half_edge;
            const verts = &self.pools.vertex;
            if (vert.halfOrNull(verts)) |half| {
                var current = half;
                while (true) {
                    if (current.leftFaceOrNull(hedges) == null) {
                        return current;
                    }
                    current = current
                        .twin(hedges)
                        .next(hedges);
                    if (current.inner == half.inner) break;
                }
            } else return null;
            return error.VertexSurronded;
        }

        /// Given a vertex finds the first inbound half edge without a left face.
        /// Returns error in the case that there is not a half edge satisfying this property
        pub fn findFirstFreeInboundHalfEdgeOrNull(self: *Self, vert: Vertex.Handle) !?HalfEdge.Handle {
            const hedges = &self.hedges;
            const verts = &self.verts;
            if (vert.halfOrNull(verts)) |half| {
                const start = half.twin(hedges);
                var current: HalfEdge.Handle = start;
                while (true) {
                    if (current.leftFaceOrNull(hedges) == null) {
                        return current;
                    }
                    current = current
                        .next(hedges)
                        .twin(hedges);
                    if (current.inner == start.inner) break;
                }
            } else return null;
            return error.VertexSurronded;
        }

        pub fn findBoundaryHalfEdgeOrNull(self: *Self, from: Vertex.Handle, to: Vertex.Handle) ?HalfEdge.Handle {
            const verts = &self.verts;
            const hedges = &self.hedges;
            if (from.halfOrNull(verts)) |from_half| {
                if (to.halfOrNull(verts)) |_| {
                    const start = from_half;
                    var current = start;
                    while (true) {
                        if (current.destination(hedges).inner == to.inner) {
                            return if (current.leftFaceOrNull(hedges) == null)
                                current
                            else
                                null;
                        }
                        current = current
                            .twin(hedges)
                            .next(hedges);
                        if (current.inner == start.inner) break;
                    }
                } else return null;
            } else return null;
            return null;
        }

        /// Simple utility to link together two half edges like links in a chain
        pub fn chainHalfEdges(self: *Self, a: HalfEdge.Handle, b: HalfEdge.Handle) void {
            const hedge_pool = &self.hedges;
            a.deref(hedge_pool).next = b;
            b.deref(hedge_pool).prev = a;
        }

        pub fn addVertex(self: *Self, allocator: Allocator, data: VertexData) !Vertex.Handle {
            const verts = &self.verts;
            const vert = try Vertex.Handle.alloc(allocator, verts);
            vert.deref(verts).* = .{
                .half = null,
                .data = data,
            };
            return vert;
        }

        /// Remove all of the edges connected to the vertex.
        /// Remove all the faces incident to those edges
        /// Deallocate the vertex.
        pub fn removeVertex(self: *Self, allocator: Allocator, vertex: Vertex.Handle) !void {
            const verts = &self.verts;
            const hedges = &self.hedges;

            if (vertex.halfOrNull(verts)) |start| {
                var to_remove: std.ArrayListUnmanaged(Edge.Handle) = .empty;
                defer to_remove.deinit(allocator);
                var current = start;
                while (true) {
                    try to_remove.append(allocator, current.parent(hedges));
                    current = current.twin(hedges).next(hedges);
                    if (current.inner == start.inner) break;
                }

                for (to_remove.items) |edge|
                    try self.removeEdge(allocator, edge);
            }

            try vertex.free(allocator, verts);
        }

        ///- An edge can be added between vertices A and B iff A and B are free.
        ///- A vertex is free if it has at least one free outgoing half-edge.
        ///- A half-edge is free if half-edge.left_face == null.
        ///- Only one edge is allowed between A and B.
        ///- A cannot equal B
        pub fn addEdge(self: *Self, allocator: Allocator, from: Vertex.Handle, to: Vertex.Handle, data: ?EdgeData) !Edge.Handle {
            const verts = &self.verts;
            const edges = &self.edges;
            const hedges = &self.hedges;

            if (to.inner == from.inner) return error.DuplicateVertices;
            if (from.halfOrNull(verts)) |half|
                if (half.destination(hedges).inner == to.inner) return error.EdgeAlreadyExists;

            const maybe_inbound_from = try self.findFirstFreeInboundHalfEdgeOrNull(from);
            const maybe_inbound_to = try self.findFirstFreeInboundHalfEdgeOrNull(to);

            // Preallocate
            const edge: Edge.Handle = try .alloc(allocator, edges);
            errdefer edge.free(allocator, edges) catch {};

            const from_half: HalfEdge.Handle = try .alloc(allocator, hedges);
            errdefer from_half.free(allocator, hedges) catch {};

            const to_half: HalfEdge.Handle = try .alloc(allocator, hedges);
            errdefer to_half.free(allocator, hedges) catch {};

            // Initialize Edge
            const edge_ptr = edge.deref(edges);
            edge_ptr.half = from_half;
            if (data) |d| edge_ptr.data = d;
            // Initialize the half edges
            const from_half_ptr = from_half.deref(hedges);
            from_half_ptr.parent = edge;
            from_half_ptr.left_face = null;
            from_half_ptr.origin = from;
            from_half_ptr.twin = to_half;
            const to_half_ptr = to_half.deref(hedges);
            to_half_ptr.parent = edge;
            to_half_ptr.left_face = null;
            to_half_ptr.origin = to;
            to_half_ptr.twin = from_half;

            if (maybe_inbound_from) |inbound_from| {
                const outbound_from = inbound_from.next(hedges);
                self.chainHalfEdges(inbound_from, from_half);
                self.chainHalfEdges(to_half, outbound_from);
            } else {
                // the vertex is isolated
                from.deref(verts).half = from_half;
                self.chainHalfEdges(to_half, from_half);
            }

            if (maybe_inbound_to) |inbound_to| {
                const outbound_to = inbound_to.next(hedges);
                self.chainHalfEdges(inbound_to, to_half);
                self.chainHalfEdges(from_half, outbound_to);
            } else {
                // the vertex is isolated
                to.deref(verts).half = to_half;
                self.chainHalfEdges(from_half, to_half);
            }

            return edge;
        }

        /// Removes an edge from the mesh, cleaning up after itself.
        pub fn removeEdge(self: *Self, allocator: Allocator, edge: Edge.Handle) !void {
            const edges = &self.edges;
            const hedges = &self.hedges;
            const verts = &self.verts;

            const half = edge.half(edges);
            const twin = half.twin(hedges);

            // Destroy adjacent faces
            if (half.leftFaceOrNull(hedges)) |face|
                try self.removeFace(allocator, face);

            if (twin.leftFaceOrNull(hedges)) |face|
                try self.removeFace(allocator, face);

            // Update the origin vertices
            { // Clean up half.origin
                const next_outbound = twin.next(hedges);
                const is_only_outbound_half_edge = next_outbound.inner == half.inner;

                half.origin(hedges).deref(verts).half =
                    if (is_only_outbound_half_edge) null else next_outbound;
            }
            { // Clean up twin.origin
                const next_outbound = half.next(hedges);
                const is_only_outbound_half_edge = next_outbound.inner == twin.inner;

                twin.origin(hedges).deref(verts).half =
                    if (is_only_outbound_half_edge) null else next_outbound;
            }
            // Destroy the half edges
            try half.free(allocator, hedges);
            try twin.free(allocator, hedges);
            // Destroy the edge
            try edge.free(allocator, edges);
        }

        /// Adds a face described by a set of vertices.
        /// Assumes that vertices are in either clockwise or counter clockwise order.
        pub fn addFace(self: *Self, allocator: Allocator, vertices: []const Vertex.Handle, data: FaceData) !Face.Handle {
            const edges = &self.edges;
            const faces = &self.faces;
            const hedges = &self.hedges;
            // Check that the data is valid and that the given half-edge loop can be turned into a polygon as given by the required conditions.
            if (vertices.len < 3) return error.LoopTooSmall;
            // TODO: Check that vertices are unique

            var construction_ring: std.ArrayListUnmanaged(HalfEdge.Handle) = try .initCapacity(allocator, vertices.len);
            defer construction_ring.deinit(allocator);
            for (vertices, 0..) |vert, idx| {
                const jdx = (idx + 1) % vertices.len;
                const next = vertices[jdx];
                if (self.findBoundaryHalfEdgeOrNull(vert, next)) |half| {
                    try construction_ring.append(allocator, half);
                } else {
                    // There is no connection between the vertices so we will make one
                    const e = try self.addEdge(allocator, vert, next, null);
                    try construction_ring.append(allocator, e.half(edges));
                }
            }
            std.debug.assert(construction_ring.items.len == vertices.len);

            const f = try Face.Handle.alloc(allocator, faces);
            f.deref(faces).* = .{
                .data = data,
                .half = construction_ring.items[0],
            };

            for (construction_ring.items) |half| {
                std.debug.assert(half.leftFaceOrNull(hedges) == null);
                half.deref(hedges).left_face = f;
            }

            return f;
        }

        /// Set the polygon of each boundary half-edge to null.
        /// Deallocate the polygon.
        pub fn removeFace(self: *Self, allocator: Allocator, face: Face.Handle) !void {
            const faces = &self.faces;
            const hedges = &self.hedges;

            const start = face.half(faces);
            var current = start;
            while (true) {
                current.deref(hedges).left_face = null;

                current = current.next(hedges);
                if (current.inner == start.inner) break;
            }

            try face.free(allocator, faces);
        }

        pub fn swapEdge(self: *Self, edge: Edge.Handle) !void {
            _ = self; // autofix
            _ = edge; // autofix
            @panic("TODO");
        }

        pub fn dissolveEdge(self: *Self, allocator: Allocator, edge: Edge.Handle) !void {
            _ = allocator; // autofix
            _ = self; // autofix
            _ = edge; // autofix
            @panic("TODO");
        }

        pub fn dissolveVertex(self: *Self, allocator: Allocator, vertex: Vertex.Handle) !void {
            _ = vertex; // autofix
            _ = allocator; // autofix
            _ = self; // autofix
            @panic("TODO");
        }

        pub fn collapseEdge(self: *Self, allocator: Allocator, edge: Edge.Handle) !void {
            _ = allocator; // autofix
            _ = self; // autofix
            _ = edge; // autofix
            @panic("TODO");
        }

        pub fn collapseFace(self: *Self, allocator: Allocator, face: Face.Handle) !void {
            _ = face; // autofix
            _ = allocator; // autofix
            _ = self; // autofix
            @panic("TODO");
        }

        pub fn subdivideEdge(self: *Self, allocator: Allocator, edge: Edge.Handle) !void {
            _ = edge; // autofix
            _ = allocator; // autofix
            _ = self; // autofix
            @panic("TODO");
        }

        pub fn subdivideFace(self: *Self, allocator: Allocator, face: Face.Handle) !void {
            _ = face; // autofix
            _ = allocator; // autofix
            _ = self; // autofix
            @panic("TODO");
        }

        pub fn addTriangle(
            self: *Self,
            allocator: Allocator,
            a: Vertex.Handle,
            b: Vertex.Handle,
            c: Vertex.Handle,
        ) !Face.Handle {
            _ = a; // autofix
            _ = b; // autofix
            _ = c; // autofix
            _ = allocator; // autofix
            _ = self; // autofix
            @panic("TODO");
        }
    };
}

const testing = std.testing;
const TestMesh = PolyMesh(void, void, void);

fn validateEdgeLoop(
    mesh: *TestMesh,
    start: TestMesh.HalfEdge.Handle,
    expect: struct {
        loop_size: ?usize,
        left_face: ?TestMesh.Face.Handle,
    },
) !void {
    const hedges = &mesh.pools.half_edge;
    var current = start;
    var count: usize = 0;
    while (true) {
        if (count >= 100) return error.LoopTooBig;
        if (expect.loop_size) |size| if (count >= size) return error.LoopTooBig;

        if (current.leftFaceOrNull(hedges)) |f1| {
            if (expect.left_face) |f2| {
                if (f1.inner != f2.inner) return error.LoopDiscontinuity;
            } else {
                return error.LoopDiscontinuity;
            }
        } else if (expect.left_face != null) return error.LoopDiscontinuity;

        current = current.next(hedges);
        count += 1;
        if (current.inner == start.inner) break;
    }
}

fn validateMesh(mesh: *TestMesh) !void {
    const hedges = &mesh.pools.half_edge;
    const edges = &mesh.pools.edge;
    const verts = &mesh.pools.vertex;
    const faces = &mesh.pools.face;
    { // Test HalfEdges
        var it = hedges.iterator();
        while (it.next()) |i| {
            const half: TestMesh.HalfEdge.Handle = .from(i);
            const parent = half.parent(hedges);

            try testing.expect(@intFromEnum(half.inner) < mesh.pools.half_edge.slot_list.items.len);

            // Make sure that there is no invalid half edge linkages
            try testing.expect(half.twin(hedges).twin(hedges).inner == half.inner);
            try testing.expect(half.next(hedges).prev(hedges).inner == half.inner);
            try testing.expect(half.prev(hedges).next(hedges).inner == half.inner);
            try testing.expect(half.prev(hedges).next(hedges).inner == half.inner);

            // Test Half edge - edge linkage
            const ph = parent.half(edges);
            try testing.expect(ph.inner == half.inner or ph.twin(hedges).inner == half.inner);
        }
    }

    { // Test Edges
        var it = edges.iterator();
        while (it.next()) |i| {
            const edge: TestMesh.Edge.Handle = .from(i);
            const half = edge.half(edges);

            try testing.expect(@intFromEnum(edge.inner) < mesh.pools.edge.slot_list.items.len);

            try testing.expect(half.origin(hedges).inner != half.destination(hedges).inner);
        }
    }

    { // Test Verts
        var it = verts.iterator();
        while (it.next()) |i| {
            const vert: TestMesh.Vertex.Handle = .from(i);
            if (vert.halfOrNull(verts)) |half| {
                try testing.expect(@intFromEnum(vert.inner) < mesh.pools.vertex.slot_list.items.len);

                try testing.expect(half.origin(hedges).inner == vert.inner);
                try testing.expect(half.twin(hedges).destination(hedges).inner == vert.inner);
            }
        }
    }

    { // Test Faces
        var it = faces.iterator();
        while (it.next()) |i| {
            const face: TestMesh.Face.Handle = .from(i);
            const half = face.half(faces);

            try testing.expect(@intFromEnum(face.inner) < mesh.pools.face.slot_list.items.len);
            try validateEdgeLoop(mesh, half, .{ .loop_size = null, .left_face = face });
        }
    }
}

test "vertex tests" {
    var mesh: TestMesh = .empty;
    defer mesh.deinit(testing.allocator);

    const v0 = try mesh.addVertex(testing.allocator, {});

    try testing.expect(v0.deref(&mesh.pools.vertex).half == null);
}

test "edge tests" {
    var mesh: TestMesh = .empty;
    defer mesh.deinit(testing.allocator);

    const edges = &mesh.pools.edge;
    const hedges = &mesh.pools.half_edge;
    const verts = &mesh.pools.vertex;

    //  v0 ---e0--- v1
    const v0 = try mesh.addVertex(testing.allocator, {});
    const v1 = try mesh.addVertex(testing.allocator, {});
    const e0 = try mesh.addEdge(testing.allocator, v0, v1, {});

    try validateMesh(&mesh);

    // did we properly link the vertices?
    const v1_saved = v1.halfOrNull(verts);
    {
        try testing.expect(v0.halfOrNull(verts) != null);
        try testing.expect(v1.halfOrNull(verts) != null);
        // does edge.half == v0.half?
        try testing.expect(e0.half(edges).inner == v0.halfOrNull(verts).?.inner);
        // does edge.half.twin == v1.half?
        try testing.expect(e0.half(edges).twin(hedges).inner == v1.halfOrNull(verts).?.inner);
    }
    // did we properly link half edges?
    {
        const start = e0.half(edges);
        var current = start;
        for (0..2) |_| {
            try testing.expect(current.leftFaceOrNull(hedges) == null);
            current = current.next(hedges);
        }
        try testing.expect(current.inner == start.inner);
    }

    //  v0 ---e0--- v1 ---e1--- v2
    const v2 = try mesh.addVertex(testing.allocator, {});
    const e1 = try mesh.addEdge(testing.allocator, v1, v2, {});

    try validateMesh(&mesh);

    // did we properly link the vertices?
    {
        // This operation should have left v1 unchanged
        try testing.expect(v1_saved.?.inner == v1.halfOrNull(verts).?.inner);
        try testing.expect(v2.halfOrNull(verts) != null);
        // does edge.half.twin == v2.half?
        try testing.expect(e1.half(edges).twin(hedges).inner == v2.halfOrNull(verts).?.inner);
    }
    // did we properly link half edges?
    {
        try testing.expect(e1.half(edges).next(hedges).inner == e1.half(edges).twin(hedges).inner);
        try testing.expect(e1.half(edges).next(hedges).next(hedges).inner == e0.half(edges).twin(hedges).inner);

        const start = e0.half(edges);
        var current = start;
        for (0..4) |_| {
            try testing.expect(current.leftFaceOrNull(hedges) == null);
            current = current.next(hedges);
        }
        try testing.expect(current.inner == start.inner);
    }

    //  v0 ---e0--- v1 ---e1--- v2
    //              |
    //              e2
    //              |
    //              v3
    const v3 = try mesh.addVertex(testing.allocator, {});
    const e2 = try mesh.addEdge(testing.allocator, v1, v3, {});

    try validateMesh(&mesh);
    // did we properly link the vertices?
    {
        // This operation should have left v1 unchanged
        try testing.expect(v1_saved.?.inner == v1.halfOrNull(verts).?.inner);
        try testing.expect(v3.halfOrNull(verts) != null);
        // does edge.half.twin == v3.half?
        try testing.expect(e2.half(edges).twin(hedges).inner == v3.halfOrNull(verts).?.inner);
    }
    // did we properly link half edges?
    {
        try testing.expect(e2.half(edges).next(hedges).inner == e2.half(edges).twin(hedges).inner);
        const sequence = .{ 0, 2, 2, 1, 1, 0 };

        const start = e0.half(edges);
        var current = start;
        inline for (0..6) |i| {
            try testing.expect(current.leftFaceOrNull(hedges) == null);
            try testing.expect(@intFromEnum(current.parent(hedges).inner) == sequence[i]);
            current = current.next(hedges);
        }
        try testing.expect(current.inner == start.inner);
    }

    try testing.expectError(error.EdgeAlreadyExists, mesh.addEdge(testing.allocator, v0, v1, {}));
}

test "correctly remove edge" {
    var mesh: TestMesh = .empty;
    defer mesh.deinit(testing.allocator);

    const edges = &mesh.pools.edge;
    _ = edges; // autofix
    const hedges = &mesh.pools.half_edge;
    _ = hedges; // autofix
    const verts = &mesh.pools.vertex;
    _ = verts; // autofix

    //  v0 ---e0--- v1 ---e1--- v2
    //              |
    //              e2
    //              |
    //              v3
    const v0 = try mesh.addVertex(testing.allocator, {});
    const v1 = try mesh.addVertex(testing.allocator, {});
    const v2 = try mesh.addVertex(testing.allocator, {});
    const v3 = try mesh.addVertex(testing.allocator, {});
    const e0 = try mesh.addEdge(testing.allocator, v0, v1, {});
    _ = e0;
    const e1 = try mesh.addEdge(testing.allocator, v1, v2, {});
    _ = e1;
    const e2 = try mesh.addEdge(testing.allocator, v1, v3, {});
    _ = e2; // autofix

}
