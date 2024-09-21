const std = @import("std");
const mem = std.mem;
const Allocator = mem.Allocator;

pub const Options = struct {
    panic_on_use_after_free: bool = true,
    panic_on_double_free: bool = true,
};

pub fn SlotPool(comptime T: type, options: Options) type {
    return struct {
        const Self = @This();
        pub const Handle = enum(u32) {
            _,
            pub fn deref(self: @This(), pool: *Self) *T {
                return pool.deref(self);
            }
        };
        pub const List = std.ArrayListUnmanaged(T);
        pub const FreeList = std.ArrayListUnmanaged(Handle);

        slot_list: List,
        free_list: FreeList,
        free_bitset: std.DynamicBitSetUnmanaged,

        pub const empty: Self = .{
            .slot_list = .empty,
            .free_list = .empty,
            .free_bitset = .{},
        };

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.slot_list.deinit(allocator);
            self.free_list.deinit(allocator);
            self.free_bitset.deinit(allocator);
        }

        pub fn initCapacity(allocator: Allocator, num: usize) Self {
            return .{
                .slot_list = List.initCapacity(allocator, num),
                .free_list = .empty,
                .free_bitset = std.DynamicBitSetUnmanaged.initEmpty(allocator, num),
            };
        }

        pub fn count(self: *Self) usize {
            return self.slot_list.items.len - self.free_list.items.len;
        }

        /// Allocates a Slot in the pool and returns the handle
        pub fn alloc(self: *Self, allocator: Allocator) !Handle {
            if (self.free_list.popOrNull()) |first_free| {
                self.free_bitset.unset(@intFromEnum(first_free));
                return first_free;
            } else {
                const len = self.slot_list.items.len;
                _ = try self.slot_list.addOne(allocator);
                try self.free_bitset.resize(allocator, self.free_bitset.bit_length + 1, false);
                return @enumFromInt(len);
            }
        }
        /// Allocates then initializes a slot with the passed value and returns the handle
        pub fn make(self: *Self, allocator: Allocator, value: T) !Handle {
            const handle = try self.alloc(allocator);
            self.set(handle, value);
            return handle;
        }

        /// Marks a slot as free
        pub fn free(self: *Self, allocator: Allocator, index: Handle) !void {
            if (options.panic_on_use_after_free and self.free_bitset.isSet(@intFromEnum(index)))
                @panic("Double Free");
            // Mark the slot as free
            self.free_bitset.set(@intFromEnum(index));
            // Append the slot to the free list
            try self.free_list.append(allocator, index);
        }

        /// Sets the value in the given slot
        pub fn set(self: *Self, index: Handle, value: T) void {
            if (options.panic_on_use_after_free and self.free_bitset.isSet(@intFromEnum(index)))
                @panic("Use After Free");
            self.slot_list.items[@intFromEnum(index)] = value;
        }

        /// Get the value in the given slot
        pub fn get(self: *Self, index: Handle) T {
            if (options.panic_on_use_after_free and self.free_bitset.isSet(@intFromEnum(index)))
                @panic("Use After Free");
            return self.slot_list.items[@intFromEnum(index)];
        }

        /// Get the value in the given slot
        pub fn deref(self: *Self, index: Handle) *T {
            if (options.panic_on_use_after_free and self.free_bitset.isSet(@intFromEnum(index)))
                @panic("Use After Free");
            return &self.slot_list.items[@intFromEnum(index)];
        }

        pub const ForwardIterator = struct {
            pool: *Self,
            idx: usize,

            pub fn next(self: *@This()) ?Handle {
                if (self.idx >= self.pool.slot_list.items.len) return null;

                if (!self.pool.free_bitset.isSet(self.idx)) {
                    defer self.idx += 1;
                    return @enumFromInt(self.idx);
                } else {
                    self.idx += 1;
                    return self.next();
                }
            }
            /// Reset the iterator to the initial index
            pub fn reset(self: *@This()) void {
                self.idx = 0;
            }
        };
        pub fn iterator(self: *Self) ForwardIterator {
            return .{
                .pool = self,
                .idx = 0,
            };
        }

        pub const ReverseIterator = struct {
            pool: *Self,
            idx: usize,

            pub fn next(self: *@This()) ?*T {
                if (self.idx < 1) return null;

                if (!self.pool.free_bitset.isSet(self.idx - 1)) {
                    defer self.idx -= 1;
                    return @enumFromInt(self.idx - 1);
                } else {
                    self.idx -= 1;
                    return self.next();
                }
            }

            /// Reset the iterator to the initial index
            pub fn reset(self: *@This()) void {
                self.idx = 0;
            }
        };
        pub fn reverseIterator(self: *Self) ReverseIterator {
            return .{
                .pool = self,
                .idx = self.slot_list.items.len,
            };
        }
    };
}
