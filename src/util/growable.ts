type TypedArrayCtor<T> = { new (length: number): T };

class GrowableArray<T extends Float64Array | Uint32Array | Uint8Array> {
  protected data: T;
  length = 0;
  private readonly ctor: TypedArrayCtor<T>;

  constructor(ctor: TypedArrayCtor<T>, initialCapacity: number) {
    this.ctor = ctor;
    this.data = new ctor(Math.max(16, initialCapacity));
  }

  push(value: number): void {
    if (this.length === this.data.length) {
      const next = new this.ctor(this.data.length * 2);
      next.set(this.data as never);
      this.data = next;
    }
    this.data[this.length++] = value;
  }

  get(index: number): number {
    return this.data[index]!;
  }

  set(index: number, value: number): void {
    this.data[index] = value;
  }

  view(): T {
    return this.data.subarray(0, this.length) as T;
  }
}

export class Float64List extends GrowableArray<Float64Array> {
  constructor(initialCapacity = 1024) {
    super(Float64Array, initialCapacity);
  }
}

export class Uint32List extends GrowableArray<Uint32Array> {
  constructor(initialCapacity = 1024) {
    super(Uint32Array, initialCapacity);
  }
}

export class Uint8List extends GrowableArray<Uint8Array> {
  constructor(initialCapacity = 1024) {
    super(Uint8Array, initialCapacity);
  }
}
