const std = @import("std");

pub fn GenericPolyMesh(comptime VertexData: type, comptime EdgeData: type, comptime FaceData: type) type {
    return struct {
        pub const Vertex = struct {
            const List = std.DoublyLinkedList(@This());
            const Pool = std.heap.MemoryPool(List.Node);
            half: ?*HalfEdge,
            data: VertexData,
        };
        pub const Edge = struct {
            const List = std.DoublyLinkedList(@This());
            const Pool = std.heap.MemoryPool(List.Node);
            half: *HalfEdge,
            data: EdgeData,
        };
        pub const Face = struct {
            const List = std.DoublyLinkedList(@This());
            const Pool = std.heap.MemoryPool(List.Node);
            half: *HalfEdge,
            data: FaceData,

            pub const VertexIterator = struct {
                current_handle: *Vertex,
                start_handle: *Vertex,
                start: bool = true,

                pub fn next(it: *@This()) ?*Vertex {
                    if (it.current_handle == it.start_handle and it.start == false) return null;
                    it.start = false;
                    const current = it.current_handle;
                    it.current_handle = current.incident_edge.?.next.origin;
                    return current;
                }

                /// Reset the iterator to the initial index
                pub fn reset(it: *@This()) void {
                    it.current_handle = it.start_handle;
                    it.start = true;
                }
            };

            pub fn vertexIter(self: *const @This()) VertexIterator {
                return .{
                    .current_handle = self.half_edge.origin,
                    .start_handle = self.half_edge.origin,
                };
            }
        };

        pub const HalfEdge = struct {
            const List = std.DoublyLinkedList(@This());
            const Pool = std.heap.MemoryPool(List.Node);
            origin: *Vertex,
            left_face: ?*Face,
            parent: *Edge,

            twin: *HalfEdge,
            next: *HalfEdge,
            prev: *HalfEdge,

            pub fn destination(self: *const @This()) *Vertex {
                return self.twin.origin;
            }
        };

        vertex_list: Vertex.List,
        edge_list: Edge.List,
        face_list: Face.List,
        half_edge_list: HalfEdge.List,

        vertex_pool: Vertex.Pool,
        edge_pool: Edge.Pool,
        face_pool: Face.Pool,
        half_edge_pool: HalfEdge.Pool,

        pub fn init(allocator: std.mem.Allocator) @This() {
            return .{
                .vertex_list = .{},
                .edge_list = .{},
                .face_list = .{},
                .half_edge_list = .{},

                .vertex_pool = .init(allocator),
                .edge_pool = .init(allocator),
                .face_pool = .init(allocator),
                .half_edge_pool = .init(allocator),
            };
        }

        pub fn deinit(self: *@This()) void {
            self.vertex_pool.deinit();
            self.edge_pool.deinit();
            self.face_pool.deinit();
            self.half_edge_pool.deinit();
        }

        pub fn allocVertex(self: *@This()) !*Vertex {
            const node = try self.vertex_pool.create();
            self.vertex_list.append(node);
            return &node.data;
        }

        pub fn createVertex(self: *@This(), vertex: Vertex) !*Vertex {
            const data = try self.allocVertex();
            data.* = vertex;
            return data;
        }

        pub fn allocEdge(self: *@This()) !*Edge {
            const node = try self.edge_pool.create();
            self.edge_list.append(node);
            return &node.data;
        }

        pub fn createEdge(self: *@This(), edge: Edge) !*Edge {
            const data = try self.allocEdge();
            data.* = edge;
            return data;
        }

        pub fn allocFace(self: *@This()) !*Face {
            const node = try self.face_pool.create();
            self.face_list.append(node);
            return &node.data;
        }

        pub fn createFace(self: *@This(), face: Face) !*Face {
            const data = try self.allocFace();
            data.* = face;
            return data;
        }

        pub fn allocHalfEdge(self: *@This()) !*HalfEdge {
            const node = try self.half_edge_pool.create();
            self.half_edge_list.append(node);
            return &node.data;
        }

        pub fn createHalfEdge(self: *@This(), half_edge: HalfEdge) !*HalfEdge {
            const data = try self.allocHalfEdge();
            data.* = half_edge;
            return data;
        }

        pub fn findFreeIncident(vertex: *Vertex) ?*HalfEdge {
            std.debug.assert(vertex.half != null);

            const begin = vertex.half.?.twin;
            var current = begin;
            while (true) {
                if (current.left_face == null) {
                    return current;
                }
                current = current.next.twin;

                if (current == begin) break;
            }

            return null;
        }

        pub fn findFreeIncidentRange(vertex: *Vertex, starting_from: *HalfEdge, and_before: *HalfEdge) ?*HalfEdge {
            std.debug.assert(vertex.half != null);
            std.debug.assert(starting_from.destination() == vertex);
            std.debug.assert(and_before.destination() == vertex);

            if (and_before == starting_from) {
                return null;
            }

            var current: *HalfEdge = starting_from;
            while (true) {
                if (current.left_face == null) {
                    return current;
                }
                current = current.next.twin;

                if (current == and_before) break;
            }

            return null;
        }

        /// Simple utility to link together two half edges like links in a chain
        pub fn chainHalfEdges(a: *HalfEdge, b: *HalfEdge) void {
            a.next = b;
            b.prev = a;
        }

        pub fn makeAdjacent(in: *HalfEdge, out: *HalfEdge) bool {
            if (in.next == out) {
                // Adjacency is already correct
                return true;
            } else if (findFreeIncidentRange(
                out.origin,
                out.twin,
                in,
            )) |g| {
                const h = g.next;
                const b = in.next;
                const d = out.prev;
                chainHalfEdges(in, out);

                // Fixup the local neighbourhood
                chainHalfEdges(g, b);
                chainHalfEdges(d, h);

                return true;
            } else {
                // There is no such half-edge.
                return false;
            }
        }

        pub fn addVertex(self: *@This(), data: VertexData) !*Vertex {
            return try self.createVertex(.{
                .half = null,
                .data = data,
            });
        }

        /// Create a free floating edge, if `to` or `from` are provided,
        /// the edge will come pre linked to them
        pub fn rawAddEdge(self: *@This(), maybe_from: ?*Vertex, maybe_to: ?*Vertex, data: EdgeData, try_link_vertex: bool) !*Edge {
            // Allocate data
            const edge = try self.allocEdge();
            const edge_from_half = try self.allocHalfEdge();
            const edge_to_half = try self.allocHalfEdge();

            // Initialize data
            edge.data = data;
            edge.half = edge_from_half;

            edge_from_half.next = edge_to_half;
            edge_from_half.prev = edge_to_half;
            edge_from_half.twin = edge_to_half;
            edge_from_half.parent = edge;
            edge_from_half.left_face = null;
            if (maybe_from) |from| {
                // link hedge to vertex
                edge_from_half.origin = from;
                // Link vertex to hedge
                if (try_link_vertex and from.half == null) {
                    from.half = edge_from_half;
                }
            }

            edge_to_half.next = edge_from_half;
            edge_to_half.prev = edge_from_half;
            edge_to_half.twin = edge_from_half;
            edge_to_half.parent = edge;
            edge_to_half.left_face = null;
            if (maybe_to) |to| {
                // link hedge to vertex
                edge_to_half.origin = to;
                // Link vertex to hedge
                if (try_link_vertex and to.half == null) {
                    to.half = edge_to_half;
                }
            }

            return edge;
        }

        ///- An edge can be added between vertices A and B iff A and B are free.
        ///- A vertex is free if it has at least one free outgoing half-edge.
        ///- A half-edge is free if half-edge.left_face == null.
        ///- Only one edge is allowed between A and B.
        ///- A cannot equal B
        pub fn addEdge(self: *@This(), from: *Vertex, to: *Vertex, data: EdgeData) !*Edge {
            std.debug.assert(to != from);

            const edge = self.rawAddEdge(from, to, data, false);
            const edge_from_half = edge.half;
            const edge_to_half = edge.half.twin;

            // Check to see if it is possible to create the edge
            var from_free_half: *HalfEdge = undefined;
            var to_free_half: *HalfEdge = undefined;
            if (from.half != null) {
                if (findFreeIncident(from)) |free_half| {
                    from_free_half = free_half;
                } else {
                    return error.ReservedVertexA;
                }
            }
            if (to.half != null) {
                if (findFreeIncident(to)) |free_half| {
                    to_free_half = free_half;
                } else {
                    return error.ReservedVertexB;
                }
            }

            // Link the from-side of the edge.

            if (from.half == null) {
                from.half = edge_from_half;
            } else {
                // store the half-edges going in and out of Vertex A
                const from_in = from_free_half;
                const from_out = from_in.next;

                chainHalfEdges(from_in, edge_from_half);
                chainHalfEdges(edge_to_half, from_out);
            }

            // Link the to-side of the edge.

            if (to.half == null) {
                to.half = edge_to_half;
            } else {
                const to_in = to_free_half;
                const to_out = to_in.next;

                chainHalfEdges(to_in, edge_to_half);
                chainHalfEdges(edge_from_half, to_out);
            }

            try self.all_the_crap.append(.{ .e = edge });
            return edge;
        }

        pub fn addFace(self: *@This(), loop: []*HalfEdge, data: FaceData) !*Face {
            // Perform Precondition Checks
            if (loop.len <= 0) return error.EmptyLoop;
            for (loop, 0..) |current, idx| {
                const next_idx = (idx + 1) % loop.len;
                const next = loop[next_idx];

                if (current.destination() != next.origin) return error.MalformedLoop;
                if (current.left_face != null) return error.NonManifoldLoop;
            }

            // Try to reorder the links to get a proper orientation
            for (loop, 0..) |current, idx| {
                const next_idx = (idx + 1) % loop.len;
                const next = loop[next_idx];

                if (!makeAdjacent(current, next)) return error.NonManifoldLoop;
            }

            // Create the face
            const face = try self.faces.create();

            face.* = .{
                .half = loop[0],
                .data = data,
            };

            // Link the half-edges to the face
            for (loop) |current| {
                current.left_face = face;
            }

            try self.all_the_crap.append(.{ .f = face });
            return face;
        }

        /// Remove all of the edges connected to the vertex.
        /// Deallocate the vertex.
        pub fn removeVertex(self: *@This(), vertex: *Vertex) void {
            if (vertex.half) |start| {
                var current: *HalfEdge = undefined;
                var next: *HalfEdge = start;
                while (true) {
                    current = next;
                    next = next.twin.next;

                    if (next.parent == current.parent) {
                        next = next.twin.next;
                    }
                    self.removeEdge(current.parent);

                    if (current == next) break;
                }
            }
            self.verts.destroy(vertex);
        }

        /// Remove all of the polygons connected to the edge.
        /// Link the half-edges of the edge off from the mesh.
        /// Deallocate the edge and its half-edges.
        pub fn removeEdge(self: *@This(), edge: *Edge) void {
            const edge_from_half = edge.half;
            const edge_to_half = edge_from_half.twin;

            if (edge_from_half.left_face) |face|
                self.removeFace(face);
            if (edge_to_half.left_face) |face|
                self.removeFace(face);

            // Link Edge From Side
            const from_vertex = edge_from_half.origin;
            const from_in = edge_from_half.prev;
            const from_out = edge_to_half.next;

            if (from_vertex.half == edge_from_half) {
                if (from_out == edge_from_half) {
                    from_vertex.half = null;
                } else {
                    from_vertex.half = from_out;
                }
            }
            chainHalfEdges(from_in, from_out);

            // Link Edge To Side
            const to_vertex = edge_to_half.origin;
            const to_in = edge_to_half.prev;
            const to_out = edge_from_half.next;

            if (to_vertex.half == edge_to_half) {
                if (to_out == edge_to_half) {
                    to_vertex.half = null;
                } else {
                    to_vertex.half = to_out;
                }
            }
            chainHalfEdges(to_in, to_out);

            self.half_edges.destroy(edge_from_half);
            self.half_edges.destroy(edge_to_half);
            self.edges.destroy(edge);
        }

        /// Set the polygon of each boundary half-edge to null.
        /// Deallocate the polygon.
        pub fn removeFace(self: *@This(), face: *Face) void {
            const begin: *HalfEdge = face.half;
            var current: *HalfEdge = begin;
            while (true) {
                current.left_face = null;
                current = current.next;

                if (current == begin) break;
            }

            self.faces.destroy(face);
        }
    };
}
