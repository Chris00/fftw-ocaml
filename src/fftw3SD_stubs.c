/* Precision independent wrappers */

#define PLAN_VAL(v) (* (FFTW(plan) *) Data_custom_val(v))

#define ALLOC_PLAN()                                            \
  alloc_custom(&FFTW(caml_plan_ops), sizeof(FFTW(plan)),        \
               sizeof(FFTW(plan)), 500 * sizeof(FFTW(plan)))


static void FFTW(caml_plan_finalize)(value v)
{
  FFTW(destroy_plan)(PLAN_VAL(v));
}

static int FFTW(caml_plan_compare)(value v1, value v2)
{
  CAMLparam2(v1, v2);
  /* fftw_plan is a pointer, just compare those of v1 and v2 */
  if(PLAN_VAL(v1) < PLAN_VAL(v2)) CAMLreturn(Val_int(-1));
  else if(PLAN_VAL(v1) > PLAN_VAL(v2)) CAMLreturn(Val_int(1));
  else CAMLreturn(Val_int(0));
}

/* compare v1 v2 = 0 ==> hash(v1) = hash(v2) */
static long FFTW(caml_plan_hash)(value plan)
{
  CAMLparam1(plan);
  /* We do not know much about the plan internals, just return the
     pointer value as hash. */
  CAMLreturn((long) PLAN_VAL(plan));
}


static struct custom_operations FFTW(caml_plan_ops) = {
  "fftw3_plan", /* identifier for serialization and deserialization */
  &(FFTW(caml_plan_finalize)),
  &(FFTW(caml_plan_compare)),
  &(FFTW(caml_plan_hash)),
  custom_serialize_default,
  custom_deserialize_default
};


/* Executing plans
 ***********************************************************************/

CAMLexport value FFTW(ocaml_execute)(value p)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute)(PLAN_VAL(p));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW(ocaml_execute_dft)(value p, value i, value o)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft)(PLAN_VAL(p), Data_bigarray_val(i), Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW(ocaml_execute_split_dft)(value p,
                                               value ri, value ii,
                                               value ro, value io)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft)(PLAN_VAL(p),
                          Data_bigarray_val(ri), Data_bigarray_val(ii),
                          Data_bigarray_val(ro), Data_bigarray_val(io));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW(ocaml_execute_dft_r2c)(value p, value i, value o)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft_r2c)(PLAN_VAL(p), Data_bigarray_val(i),
                        Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW(ocaml_execute_split_dft_r2c)(value p,
                                                   value i,
                                                   value ro, value io)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft_r2c)(PLAN_VAL(p),
                              Data_bigarray_val(i),
                              Data_bigarray_val(ro), Data_bigarray_val(io));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW(ocaml_execute_dft_c2r)(value p, value i, value o)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_dft_c2r)(PLAN_VAL(p), Data_bigarray_val(i),
                        Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}

CAMLexport value FFTW(ocaml_execute_split_dft_c2r)(value p,
                                                   value ri, value ii,
                                                   value o)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_split_dft_c2r)(PLAN_VAL(p),
                              Data_bigarray_val(ri), Data_bigarray_val(ii),
                              Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


CAMLexport value FFTW(ocaml_execute_r2r)(value p, value i, value o)
{
  /* noalloc */
  enter_blocking_section();  /* Allow other threads */
  FFTW(execute_r2r)(PLAN_VAL(p), Data_bigarray_val(i), Data_bigarray_val(o));
  leave_blocking_section();  /* Disallow other threads */
  return(Val_unit);
}


/* Normalizing transforms
 ***********************************************************************/

CAMLexport
value FFTW(ocaml_normalize)(value va, /* array */
                            value vofs, value vstride, value vdim,
                            value vfactor)
{
  /* noalloc */
  FFTW(complex) *a = Data_bigarray_val(va) + Int_val(vofs);
  
  enter_blocking_section();  /* Allow other threads */
/* FIXME: code missing */
  leave_blocking_section();  /* Disallow other threads */
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
value FFTW(ocaml_guru_dft)(value vi, value vo,  value sign, value flags,
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
value FFTW(ocaml_guru_dft_bc)(value * argv, int argn)
{
  return FFTW(ocaml_guru_dft)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10], argv[11]);
}


CAMLexport
value FFTW(ocaml_guru_r2c)(value vi, value vo, value flags,
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
value FFTW(ocaml_guru_r2c_bc)(value * argv, int argn)
{
  return FFTW(ocaml_guru_r2c)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10]);
}



CAMLexport
value FFTW(ocaml_guru_c2r)(value vi, value vo, value flags,
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
value FFTW(ocaml_guru_c2r_bc)(value * argv, int argn)
{
  return FFTW(ocaml_guru_c2r)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10]);
}


CAMLexport
value FFTW(ocaml_guru_r2r)(value vi, value vo, value vkind, value flags,
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

  for(k = 0; k < rank; k++)
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
value FFTW(ocaml_guru_r2r_bc)(value * argv, int argn)
{
  return FFTW(ocaml_guru_r2r)(argv[0], argv[1], argv[2], argv[3], argv[4],
                              argv[5], argv[6], argv[7], argv[8], argv[9],
                              argv[10], argv[11]);
}
