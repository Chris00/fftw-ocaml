/* Verbatim copy of functions in "bigarray_stubs.c" of the OCaml
   distribution that are needed for the specially allocated fftw
   arrays (not defined here).  Everything must be declared
   "static". */

/* From ocaml/byterun/compare.h */
CAMLextern int caml_compare_unordered;


static uintnat caml_ba_num_elts(struct caml_bigarray * b)
{
  uintnat num_elts;
  int i;
  num_elts = 1;
  for (i = 0; i < b->num_dims; i++) num_elts = num_elts * b->dim[i];
  return num_elts;
}

static int caml_ba_element_size[] =
{ 4 /*FLOAT32*/, 8 /*FLOAT64*/,
  1 /*SINT8*/, 1 /*UINT8*/,
  2 /*SINT16*/, 2 /*UINT16*/,
  4 /*INT32*/, 8 /*INT64*/,
  sizeof(value) /*CAML_INT*/, sizeof(value) /*NATIVE_INT*/,
  8 /*COMPLEX32*/, 16 /*COMPLEX64*/
};

static uintnat
caml_ba_multov(uintnat a, uintnat b, int * overflow)
{
#define HALF_SIZE (sizeof(uintnat) * 4)
#define HALF_MASK (((uintnat)1 << HALF_SIZE) - 1)
#define LOW_HALF(x) ((x) & HALF_MASK)
#define HIGH_HALF(x) ((x) >> HALF_SIZE)
  /* Cut in half words */
  uintnat al = LOW_HALF(a);
  uintnat ah = HIGH_HALF(a);
  uintnat bl = LOW_HALF(b);
  uintnat bh = HIGH_HALF(b);
  /* Exact product is:
              al * bl
           +  ah * bl  << HALF_SIZE
           +  al * bh  << HALF_SIZE
           +  ah * bh  << 2*HALF_SIZE
     Overflow occurs if:
        ah * bh is not 0, i.e. ah != 0 and bh != 0
     OR ah * bl has high half != 0
     OR ah * bl has high half != 0
     OR the sum al * bl + LOW_HALF(ah * bl) << HALF_SIZE
                        + LOW_HALF(al * bh) << HALF_SIZE overflows.
     This sum is equal to p = (a * b) modulo word size. */
  uintnat p1 = al * bh;
  uintnat p2 = ah * bl;
  uintnat p = a * b;
  if (ah != 0 && bh != 0) *overflow = 1;
  if (HIGH_HALF(p1) != 0 || HIGH_HALF(p2) != 0) *overflow = 1;
  p1 <<= HALF_SIZE;
  p2 <<= HALF_SIZE;
  p1 += p2;
  if (p < p1 || p1 < p2) *overflow = 1; /* overflow in sums */
  return p;
#undef HALF_SIZE
#undef LOW_HALF
#undef HIGH_HALF
}


#define CAML_BA_MAX_MEMORY 256*1024*1024


static int caml_ba_compare(value v1, value v2)
{
  struct caml_bigarray * b1 = Bigarray_val(v1);
  struct caml_bigarray * b2 = Bigarray_val(v2);
  uintnat n, num_elts;
  int i;

  /* Compare number of dimensions */
  if (b1->num_dims != b2->num_dims) return b2->num_dims - b1->num_dims;
  /* Same number of dimensions: compare dimensions lexicographically */
  for (i = 0; i < b1->num_dims; i++) {
    intnat d1 = b1->dim[i];
    intnat d2 = b2->dim[i];
    if (d1 != d2) return d1 < d2 ? -1 : 1;
  }
  /* Same dimensions: compare contents lexicographically */
  num_elts = caml_ba_num_elts(b1);

#define DO_INTEGER_COMPARISON(type) \
  { type * p1 = b1->data; type * p2 = b2->data; \
    for (n = 0; n < num_elts; n++) { \
      type e1 = *p1++; type e2 = *p2++; \
      if (e1 < e2) return -1; \
      if (e1 > e2) return 1; \
    } \
    return 0; \
  }
#define DO_FLOAT_COMPARISON(type) \
  { type * p1 = b1->data; type * p2 = b2->data; \
    for (n = 0; n < num_elts; n++) { \
      type e1 = *p1++; type e2 = *p2++; \
      if (e1 < e2) return -1; \
      if (e1 > e2) return 1; \
      if (e1 != e2) { \
        compare_unordered = 1; \
        if (e1 == e1) return 1; \
        if (e2 == e2) return -1; \
      } \
    } \
    return 0; \
  }

  switch (b1->flags & BIGARRAY_KIND_MASK) {
  case BIGARRAY_COMPLEX32:
    num_elts *= 2; /*fallthrough*/
  case BIGARRAY_FLOAT32:
    DO_FLOAT_COMPARISON(float);
  case BIGARRAY_COMPLEX64:
    num_elts *= 2; /*fallthrough*/
  case BIGARRAY_FLOAT64:
    DO_FLOAT_COMPARISON(double);
  case BIGARRAY_SINT8:
    DO_INTEGER_COMPARISON(int8);
  case BIGARRAY_UINT8:
    DO_INTEGER_COMPARISON(uint8);
  case BIGARRAY_SINT16:
    DO_INTEGER_COMPARISON(int16);
  case BIGARRAY_UINT16:
    DO_INTEGER_COMPARISON(uint16);
  case BIGARRAY_INT32:
    DO_INTEGER_COMPARISON(int32);
  case BIGARRAY_INT64:
#ifdef ARCH_INT64_TYPE
    DO_INTEGER_COMPARISON(int64);
#else
    { int64 * p1 = b1->data; int64 * p2 = b2->data;
      for (n = 0; n < num_elts; n++) {
        int64 e1 = *p1++; int64 e2 = *p2++;
        if ((int32)e1.h > (int32)e2.h) return 1;
        if ((int32)e1.h < (int32)e2.h) return -1;
        if (e1.l > e2.l) return 1;
        if (e1.l < e2.l) return -1;
      }
      return 0;
    }
#endif
  case BIGARRAY_CAML_INT:
  case BIGARRAY_NATIVE_INT:
    DO_INTEGER_COMPARISON(intnat);
  default:
    Assert(0);
    return 0;                   /* should not happen */
  }
#undef DO_INTEGER_COMPARISON
#undef DO_FLOAT_COMPARISON
}

/* Hashing of a bigarray */

static intnat caml_ba_hash(value v)
{
  struct caml_bigarray * b = Bigarray_val(v);
  intnat num_elts, n, h;
  int i;

  num_elts = 1;
  for (i = 0; i < b->num_dims; i++) num_elts = num_elts * b->dim[i];
  if (num_elts >= 50) num_elts = 50;
  h = 0;

#define COMBINE(h,v) ((h << 4) + h + (v))

  switch (b->flags & BIGARRAY_KIND_MASK) {
  case BIGARRAY_SINT8:
  case BIGARRAY_UINT8: {
    uint8 * p = b->data;
    for (n = 0; n < num_elts; n++) h = COMBINE(h, *p++);
    break;
  }
  case BIGARRAY_SINT16:
  case BIGARRAY_UINT16: {
    uint16 * p = b->data;
    for (n = 0; n < num_elts; n++) h = COMBINE(h, *p++);
    break;
  }
  case BIGARRAY_FLOAT32:
  case BIGARRAY_COMPLEX32:
  case BIGARRAY_INT32:
#ifndef ARCH_SIXTYFOUR
  case BIGARRAY_CAML_INT:
  case BIGARRAY_NATIVE_INT:
#endif
  {
    uint32 * p = b->data;
    for (n = 0; n < num_elts; n++) h = COMBINE(h, *p++);
    break;
  }
  case BIGARRAY_FLOAT64:
  case BIGARRAY_COMPLEX64:
  case BIGARRAY_INT64:
#ifdef ARCH_SIXTYFOUR
  case BIGARRAY_CAML_INT:
  case BIGARRAY_NATIVE_INT:
#endif
#ifdef ARCH_SIXTYFOUR
  {
    uintnat * p = b->data;
    for (n = 0; n < num_elts; n++) h = COMBINE(h, *p++);
    break;
  }
#else
  {
    uint32 * p = b->data;
    for (n = 0; n < num_elts; n++) {
#ifdef ARCH_BIG_ENDIAN
      h = COMBINE(h, p[1]); h = COMBINE(h, p[0]); p += 2;
#else
      h = COMBINE(h, p[0]); h = COMBINE(h, p[1]); p += 2;
#endif
    }
    break;
  }
#endif
  }
#undef COMBINE
  return h;
}

static void caml_ba_serialize_longarray(void * data,
                                        intnat num_elts,
                                        intnat min_val, intnat max_val)
{
#ifdef ARCH_SIXTYFOUR
  int overflow_32 = 0;
  intnat * p, n;
  for (n = 0, p = data; n < num_elts; n++, p++) {
    if (*p < min_val || *p > max_val) { overflow_32 = 1; break; }
  }
  if (overflow_32) {
    serialize_int_1(1);
    serialize_block_8(data, num_elts);
  } else {
    serialize_int_1(0);
    for (n = 0, p = data; n < num_elts; n++, p++) serialize_int_4((int32) *p);
  }
#else
  serialize_int_1(0);
  serialize_block_4(data, num_elts);
#endif
}

static void caml_ba_serialize(value v,
                              uintnat * wsize_32,
                              uintnat * wsize_64)
{
  struct caml_bigarray * b = Bigarray_val(v);
  intnat num_elts;
  int i;

  /* Serialize header information */
  serialize_int_4(b->num_dims);
  serialize_int_4(b->flags & (BIGARRAY_KIND_MASK | BIGARRAY_LAYOUT_MASK));
  for (i = 0; i < b->num_dims; i++) serialize_int_4(b->dim[i]);
  /* Compute total number of elements */
  num_elts = 1;
  for (i = 0; i < b->num_dims; i++) num_elts = num_elts * b->dim[i];
  /* Serialize elements */
  switch (b->flags & BIGARRAY_KIND_MASK) {
  case BIGARRAY_SINT8:
  case BIGARRAY_UINT8:
    serialize_block_1(b->data, num_elts); break;
  case BIGARRAY_SINT16:
  case BIGARRAY_UINT16:
    serialize_block_2(b->data, num_elts); break;
  case BIGARRAY_FLOAT32:
  case BIGARRAY_INT32:
    serialize_block_4(b->data, num_elts); break;
  case BIGARRAY_COMPLEX32:
    serialize_block_4(b->data, num_elts * 2); break;
  case BIGARRAY_FLOAT64:
  case BIGARRAY_INT64:
    serialize_block_8(b->data, num_elts); break;
  case BIGARRAY_COMPLEX64:
    serialize_block_8(b->data, num_elts * 2); break;
  case BIGARRAY_CAML_INT:
    caml_ba_serialize_longarray(b->data, num_elts, -0x40000000, 0x3FFFFFFF);
    break;
  case BIGARRAY_NATIVE_INT:
    caml_ba_serialize_longarray(b->data, num_elts, -0x80000000, 0x7FFFFFFF);
    break;
  }
  /* Compute required size in Caml heap.  Assumes struct caml_bigarray
     is exactly 4 + num_dims words */
  Assert(sizeof(struct caml_bigarray) == 5 * sizeof(value));
  *wsize_32 = (4 + b->num_dims) * 4;
  *wsize_64 = (4 + b->num_dims) * 8;
}


static void caml_ba_deserialize_longarray(void * dest, intnat num_elts)
{
  int sixty = deserialize_uint_1();
#ifdef ARCH_SIXTYFOUR
  if (sixty) {
    deserialize_block_8(dest, num_elts);
  } else {
    intnat * p, n;
    for (n = 0, p = dest; n < num_elts; n++, p++) *p = deserialize_sint_4();
  }
#else
  if (sixty)
    deserialize_error("input_value: cannot read bigarray "
                      "with 64-bit Caml ints");
  deserialize_block_4(dest, num_elts);
#endif
}
