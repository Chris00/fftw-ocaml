/* Precision independent wrappers */

#define PLAN_VAL(v) (* (FFTW(plan) *) Data_custom_val(v))

#define ALLOC_PLAN()                                            \
  alloc_custom(&FFTW_OCAML(plan_ops), sizeof(FFTW(plan)),       \
               sizeof(FFTW(plan)), 500 * sizeof(FFTW(plan)))


static void FFTW_OCAML(plan_finalize)(value v)
{
  FFTW(destroy_plan)(PLAN_VAL(v));
}

static int FFTW_OCAML(plan_compare)(value v1, value v2)
{
  CAMLparam2(v1, v2);
  /* fftw_plan is a pointer, just compare those of v1 and v2 */
  if(PLAN_VAL(v1) < PLAN_VAL(v2)) CAMLreturn(Val_int(-1));
  else if(PLAN_VAL(v1) > PLAN_VAL(v2)) CAMLreturn(Val_int(1));
  else CAMLreturn(Val_int(0));
}

/* compare v1 v2 = 0 ==> hash(v1) = hash(v2) */
static long FFTW_OCAML(plan_hash)(value plan)
{
  CAMLparam1(plan);
  /* We do not know much about the plan internals, just return the
     pointer value as hash. */
  CAMLreturn((long) PLAN_VAL(plan));
}


static struct custom_operations FFTW_OCAML(plan_ops) = {
  "fftw3_plan", /* identifier for serialization and deserialization */
  &(FFTW_OCAML(plan_finalize)),
  &(FFTW_OCAML(plan_compare)),
  &(FFTW_OCAML(plan_hash)),
  custom_serialize_default,
  custom_deserialize_default
};


/* Executing plans
 ***********************************************************************/

CAMLexport value FFTW_OCAML(execute)(value p)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute)(PLAN_VAL(p));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW_OCAML(execute_dft)(value p, value i, value o)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft)(PLAN_VAL(p), Data_bigarray_val(i),
                    Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW_OCAML(execute_split_dft)(value p,
                                               value ri, value ii,
                                               value ro, value io)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft)(PLAN_VAL(p),
                          Data_bigarray_val(ri), Data_bigarray_val(ii),
                          Data_bigarray_val(ro), Data_bigarray_val(io));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW_OCAML(execute_dft_r2c)(value p, value i, value o)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft_r2c)(PLAN_VAL(p), Data_bigarray_val(i),
                        Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW_OCAML(execute_split_dft_r2c)(value p,
                                                   value i,
                                                   value ro, value io)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft_r2c)(PLAN_VAL(p),
                              Data_bigarray_val(i),
                              Data_bigarray_val(ro), Data_bigarray_val(io));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW_OCAML(execute_dft_c2r)(value p, value i, value o)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft_c2r)(PLAN_VAL(p), Data_bigarray_val(i),
                        Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW_OCAML(execute_split_dft_c2r)(value p,
                                                   value ri, value ii,
                                                   value o)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft_c2r)(PLAN_VAL(p),
                              Data_bigarray_val(ri), Data_bigarray_val(ii),
                              Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW_OCAML(execute_r2r)(value p, value i, value o)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_r2r)(PLAN_VAL(p), Data_bigarray_val(i), Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


/* Normalizing transforms
 ***********************************************************************/

#define NORMALIZE(name, ty)                                             \
  static void normalize_ ## name(value va,                              \
                                 int ofs, value vstride, value vdim,    \
                                 ty factor)                             \
  {                                                                     \
    /* noalloc */                                                       \
    ty *a = Data_bigarray_val(va) + ofs;                                \
    ty temp;                                                            \
    int rank = Wosize_val(vstride);                                     \
    int num_dims = Bigarray_val(v)->num_dims;                           \
    int stride[MAX_NUM_DIMS], dim[MAX_NUM_DIMS];                        \
    int i, j, k, e0;                                                    \
                                                                        \
    for(k = 0; k < ; k++) {                                             \
      stride[k] = Int_val(Field(vstride, k));                           \
      dim[k] = Int_val(Field(vdim, k));                                 \
    }                                                                   \
    enter_blocking_section();  /* Allow other threads */                \
    switch(num_dims) {                                                  \
      /* FIXME: some dimensions my be used for "howmany" */             \
    case 1:                                                             \
      e0 = dim[0] * stride[0];                                          \
      for(i=0; i < e0; i += stride[0]) {                                \
        MUL_BY(a[i], factor);                                           \
      }                                                                 \
      break;                                                            \
    case 2:                                                             \
                                                                        \
      break;                                                            \
    default:                                                            \
      /* FIXME: code missing */                                         \
      caml_failwith("Fftw3." PREC ": normalize not yet implemented");   \
    }                                                                   \
    leave_blocking_section();  /* Disallow other threads */             \
  }

#define MUL_BY(x, a) x *= a
NORMALIZE(float, FLOAT)
#define ASSIGN_MUL(x, a)            \
  temp = x[0] * a[0] - x[1] * a[1]; \
  x[1] = x[0] * a[1] + x[1] * a[0]; \
  x[0] = temp
NORMALIZE(complex, FFTW(complex))

CAMLexport
value FFTW_OCAML(normalize)(value va, /* array */
                            value vofs, value vstride, value vdim,
                            value vfactor)
{
  /* noalloc */
  FFTW_RAISE_NO_FFTWF;
  switch (Bigarray_val(va)->flags & BIGARRAY_KIND_MASK) {
  case BIGARRAY_FLOAT64:
  case BIGARRAY_FLOAT32:
    normalize_float(va, Int_val(vofs), vstride, vdim, Double_val(vfactor));
    break;
  case BIGARRAY_COMPLEX64:
  case BIGARRAY_COMPLEX32:
    normalize_complex(va, Int_val(vofs), vstride, vdim, vfactor);
    break;
  default:
    caml_failwith("Fftw3." PREC ".normalize: wrong kind of bigarray. "
                  "Please report this bug.");
  }
  return(Val_unit);  
}



/* Creating ND plans
 ***********************************************************************/

/* Copy the Caml arrays to the structure required by the guru interface.
 * Single transforms may be specified with [howmany_rank = 0]. */
#define MAKE_DIMS()                                             \
  const int rank = Wosize_val(vn);                              \
  const int howmany_rank = Wosize_val(vhowmany);                \
  FFTW(iodim) dims[MAX_NUM_DIMS], howmany_dims[MAX_NUM_DIMS];   \
  int k;                                                        \
                                                                \
  for(k = 0; k < rank; k++) {                                   \
    dims[k].n = Int_val(Field(vn, k));                          \
    dims[k].is = Int_val(Field(vistride, k));                   \
    dims[k].os = Int_val(Field(vostride, k));                   \
  }                                                             \
  for(k = 0; k < howmany_rank; k++) {                           \
    howmany_dims[k].n = Int_val(Field(vhowmany, k));            \
    howmany_dims[k].is = Int_val(Field(vhowmanyi, k));          \
    howmany_dims[k].os = Int_val(Field(vhowmanyo, k));          \
  }

CAMLexport
value FFTW_OCAML(guru_dft)(value vi, value vo,  value sign, value flags,
                           value vofsi, value vofso,
                           value vn, value vistride, value vostride,
                           value vhowmany, value vhowmanyi, value vhowmanyo)
{
  CAMLparam5(vi, vo, sign, flags, vofsi);
  CAMLxparam5(vofso, vn, vistride, vostride, vhowmany);
  CAMLxparam2(vhowmanyi, vhowmanyo);
  CAMLlocal1(plan);
  FFTW(plan) p;
  FFTW(complex) *i = Data_bigarray_val(vi);
  FFTW(complex) *o = Data_bigarray_val(vo);
  MAKE_DIMS();
  
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  p = FFTW(plan_guru_dft)(rank, dims,
                          howmany_rank, howmany_dims,
                          i + Int_val(vofsi),
                          o + Int_val(vofso),
                          Int_val(sign), Int_val(flags));
  leave_blocking_section();  /* Disallow other threads */

  if (p == NULL) caml_failwith("Fftw3." PREC ".Genarray.dft");
  plan = ALLOC_PLAN();
  PLAN_VAL(plan) = p;
  CAMLreturn(plan);
}

CAMLexport
value FFTW_OCAML(guru_dft_bc)(value * argv, int argn)
{
  return FFTW_OCAML(guru_dft)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10], argv[11]);
}


CAMLexport
value FFTW_OCAML(guru_r2c)(value vi, value vo, value flags,
                           value vofsi, value vofso,
                           value vn, value vistride, value vostride,
                           value vhowmany, value vhowmanyi, value vhowmanyo)
{
  CAMLparam5(vi, vo, flags, vofsi, vofso);
  CAMLxparam5(vn, vistride, vostride, vhowmany, vhowmanyi);
  CAMLxparam1(vhowmanyo);
  CAMLlocal1(plan);
  FFTW(plan) p;
  FLOAT *i = Data_bigarray_val(vi);
  FFTW(complex) *o = Data_bigarray_val(vo);
  MAKE_DIMS();
  
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  p = FFTW(plan_guru_dft_r2c)(rank, dims,
                              howmany_rank, howmany_dims,
                              i + Int_val(vofsi),
                              o + Int_val(vofso),
                              Int_val(flags));
  leave_blocking_section();  /* Disallow other threads */

  if (p == NULL) caml_failwith("Fftw3." PREC ".Genarray.r2c");
  plan = ALLOC_PLAN();
  PLAN_VAL(plan) = p;
  CAMLreturn(plan);
}

CAMLexport
value FFTW_OCAML(guru_r2c_bc)(value * argv, int argn)
{
  return FFTW_OCAML(guru_r2c)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10]);
}



CAMLexport
value FFTW_OCAML(guru_c2r)(value vi, value vo, value flags,
                           value vofsi, value vofso,
                           value vn, value vistride, value vostride,
                           value vhowmany, value vhowmanyi, value vhowmanyo)
{
  CAMLparam5(vi, vo, flags, vofsi, vofso);
  CAMLxparam5(vn, vistride, vostride, vhowmany, vhowmanyi);
  CAMLxparam1(vhowmanyo);
  CAMLlocal1(plan);
  FFTW(plan) p;
  FFTW(complex) *i = Data_bigarray_val(vi);
  FLOAT *o = Data_bigarray_val(vo);
  MAKE_DIMS();
  
  FFTW_RAISE_NO_FFTWF;
  enter_blocking_section();  /* Allow other threads */
  p = FFTW(plan_guru_dft_c2r)(rank, dims,
                              howmany_rank, howmany_dims,
                              i + Int_val(vofsi),
                              o + Int_val(vofso),
                              Int_val(flags));
  leave_blocking_section();  /* Disallow other threads */

  if (p == NULL) caml_failwith("Fftw3." PREC ".Genarray.c2r");
  plan = ALLOC_PLAN();
  PLAN_VAL(plan) = p;
  CAMLreturn(plan);
}

CAMLexport
value FFTW_OCAML(guru_c2r_bc)(value * argv, int argn)
{
  return FFTW_OCAML(guru_c2r)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10]);
}


CAMLexport
value FFTW_OCAML(guru_r2r)(value vi, value vo, value vkind, value flags,
                           value vofsi, value vofso,
                           value vn, value vistride, value vostride,
                           value vhowmany, value vhowmanyi, value vhowmanyo)
{
  CAMLparam5(vi, vo, vkind, flags, vofsi);
  CAMLxparam5(vofso, vn, vistride, vostride, vhowmany);
  CAMLxparam2(vhowmanyi, vhowmanyo);
  CAMLlocal1(plan);
  FFTW(plan) p;
  FLOAT *i = Data_bigarray_val(vi);
  FLOAT *o = Data_bigarray_val(vo);
  FFTW(r2r_kind) kind[MAX_NUM_DIMS];
  MAKE_DIMS();

  FFTW_RAISE_NO_FFTWF;
  for(k = 0; k < rank; k++)
    /* OK because the order of "type r2r_kind" in fftw3SD.ml is in
     * sync with fftw3.h (see configure.ac). */
    kind[k] = Int_val(Field(vkind, k));
    
  enter_blocking_section();  /* Allow other threads */
  p = FFTW(plan_guru_r2r)(rank, dims,
                          howmany_rank, howmany_dims,
                          i + Int_val(vofsi),
                          o + Int_val(vofso),
                          kind, Int_val(flags));
  leave_blocking_section();  /* Disallow other threads */

  if (p == NULL) caml_failwith("Fftw3." PREC ".Genarray.r2r");
  plan = ALLOC_PLAN();
  PLAN_VAL(plan) = p;
  CAMLreturn(plan);
}

CAMLexport
value FFTW_OCAML(guru_r2r_bc)(value * argv, int argn)
{
  return FFTW_OCAML(guru_r2r)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10], argv[11]);
}
