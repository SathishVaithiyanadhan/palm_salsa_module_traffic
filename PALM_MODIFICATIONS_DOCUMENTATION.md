# PALM Modifications Documentation

> **Folder**: `JOBS/salsa_tra/USER_CODE/`
> **Date**: 2026-08-13
> **Author**: Sathish Kumar Vaithiyanadhan
> **Purpose**: Fix floating overflow crashes + Rosenbrock integrator failure + restart MPI deadlock (non-convergence)

---

## Quick Reference: Crash → Fix Mapping

| Crash / Issue | Backtrace | Fix | File Changed |
|---|---|---|---|
| `forrtl: error (72): floating overflow` | `data_output_mask.f90:1038` | NaN→0.0 before NF90_PUT_VAR | `data_output_mask.f90` |
| `forrtl: error (72): floating overflow` | `data_output_3d.f90:796` | NaN→0.0 before NF90_PUT_VAR | `data_output_3d.f90` |
| `forrtl: error (72): floating overflow` | Any `data_output_*.f90` | NaN→0.0 (all output types) | All 5 output files |
| `Step size too small` (Rosenbrock) | Rosenbrock ODE solver | 1.0 molec/cm³ floor | `chemistry_model_mod.f90` |
| Restart array overflow | `restart_data_mpi_io_mod` | `max_nr_arrays`: 600 → 5000 | `restart_data_mpi_io_mod.f90` |
| Restart "hang" (MPI deadlock) | `equilibration` / `distr_update` | non-fatal break-out + `SALSA-WARN` | `salsa_mod.f90` |
| `forrtl: error (65): floating invalid` | `salsa_mod.f90:18046` (`set_nest_mass`) | NaN/Inf/non-positive operand guard → 0.0 + log | `salsa_mod.f90` |

---

## File 1: `data_output_mask.f90`

**Purpose**: Writes masked output variables (`data_output_masks` in namelist)

### Change 1.1 — Add `nan_count` variable

**Line**: ~178 (after `INTEGER(iwp) :: n`)
**Why**: Counter for how many grid cells have NaN/Inf values

**Original** (source):                     | **Modified** (USER_CODE):
---                                       | ---
```fortran                                 | ```fortran
INTEGER(iwp) ::  n                        | INTEGER(iwp) ::  n
#if defined( __parallel )                 | INTEGER(iwp) ::  nan_count
INTEGER(iwp) ::  ind(6)                   | #if defined( __parallel )
```                                        | INTEGER(iwp) ::  ind(6)
```                                        | ```

### Change 1.2 — NaN protection (Parallel netCDF4 path, `local_pf`)

**Lines**: ~986-1001
**Array**: `local_pf` — local masked data on each PE
**When**: Any cell has NaN, +Inf, or -Inf

**Original** (~line 988):
```fortran
          THEN
             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                                     &
                                     id_var_domask(mid,av,ivar), local_pf,                    &
```

**Modified** (~line 986):
```fortran
          THEN
!
!--          Replace NaN/Inf with 0.0 to prevent NetCDF write failure (floating overflow).
             nan_count = COUNT( local_pf /= local_pf ) + COUNT( local_pf > 1.0E30_wp ) +     &
                          COUNT( local_pf < -1.0E30_wp )
             IF ( nan_count > 0 )  THEN
                WRITE( message_string, * ) 'data_output_mask: NaN/Inf detected in variable ''',&
                                           TRIM( domask(mid,av,ivar) ), ''' at timestep ',     &
                                           domask_time_count(mid,av), ', mask ', mid,         &
                                           ' (', nan_count,                                   &
                                           ' cells), replaced with 0.0 for masked output'
                CALL message( 'data_output_mask', 'PAC1000', 0, 1, 0, 6, 0 )
                WHERE( local_pf /= local_pf .OR. local_pf > 1.0E30_wp .OR.                   &
                       local_pf < -1.0E30_wp )                                                &
                   local_pf = 0.0_wp
             ENDIF
             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                                     &
                                     id_var_domask(mid,av,ivar), local_pf,                    &
```

### Change 1.3 — NaN protection (Parallel conventional path, `total_pf`) ← **CRASH SITE**

**Lines**: ~1056-1071
**Array**: `total_pf` — gathered masked data on PE0
**Why**: This is the PRIMARY crash site — `data_output_mask.f90:1038` in backtrace

**Original** (~line 1045):
```fortran
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                                     &
                                     id_var_domask(mid,av,ivar), total_pf,                    &
```

**Modified** (~line 1056):
```fortran
             ENDDO

!
!--          Replace NaN/Inf with 0.0 to prevent NetCDF write failure (floating overflow).
             nan_count = COUNT( total_pf /= total_pf ) + COUNT( total_pf > 1.0E30_wp ) +     &
                          COUNT( total_pf < -1.0E30_wp )
             IF ( nan_count > 0 )  THEN
                WRITE( message_string, * ) 'data_output_mask: NaN/Inf detected in variable ''',&
                                           TRIM( domask(mid,av,ivar) ), ''' at timestep ',     &
                                           domask_time_count(mid,av), ', mask ', mid,         &
                                           ' (', nan_count,                                   &
                                           ' cells), replaced with 0.0 for masked output'
                CALL message( 'data_output_mask', 'PAC1000', 0, 1, 0, 6, 0 )
                WHERE( total_pf /= total_pf .OR. total_pf > 1.0E30_wp .OR.                    &
                       total_pf < -1.0E30_wp )                                                 &
                   total_pf = 0.0_wp
             ENDIF
             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                                     &
                                     id_var_domask(mid,av,ivar), total_pf,                    &
```

### Change 1.4 — NaN protection (Serial path, `local_pf`)

**Lines**: ~1110-1124
**Array**: `local_pf` — serial (single processor) execution

Same pattern as 1.2 in the `#else` (non-parallel) branch.

---

## File 2: `data_output_3d.f90`

**Purpose**: Writes 3D volume output (`data_output` in namelist)

### Change 2.1 — Add `nan_count` variable

**Line**: ~197 (after `INTEGER(iwp) :: nzt_do`)

### Change 2.2 — NaN protection (Parallel netCDF4, `tr_out`)

**Lines**: ~834-847 | **Array**: `tr_out` — transposed 3D array

**Original**:                                    | **Modified**:
---                                              | ---
```fortran                                        | ```fortran
!$ACC UPDATE HOST(tr_out)                        | !$ACC UPDATE HOST(tr_out)
nc_stat = NF90_PUT_VAR(... tr_out ...)            | nan_count = COUNT( tr_out /= tr_out ) + &
```                                               |              COUNT( tr_out > 1.0E30_wp ) + &
                                                  |              COUNT( tr_out < -1.0E30_wp )
                                                  | IF ( nan_count > 0 )  THEN
                                                  |    WRITE( message_string, * ) 'data_output_3d: ...'
                                                  |    CALL message( ... )
                                                  |    WHERE( tr_out /= ... ) tr_out = 0.0_wp
                                                  | ENDIF
                                                  | nc_stat = NF90_PUT_VAR(... tr_out ...)
                                                  | ```

### Change 2.3 — NaN protection (Parallel conventional, `local_pf`)

**Lines**: ~870-882 | Same pattern on `local_pf`

### Change 2.4 — NaN protection (Serial, `local_pf`)

**Lines**: ~900-910 | Same pattern in `#else` (non-parallel) branch

---

## File 3: `data_output_2d.f90`

**Purpose**: Writes 2D cross-sections (`_xy`, `_xz`, `_yz` suffixes)

### Change 3.1 — Add `nan_count` variable

**Line**: ~193 (after `INTEGER(iwp) :: s_ind`)

### Changes 3.2–3.10 — NaN protection (9 locations)

| # | Lines | Array Protected | Slice | Path |
|---|---|---|---|---|
| 1 | ~1626-1638 | `total_2d` | XY | Parallel (PE0 gather) |
| 2 | ~1677-1690 | `local_2d` | XY | Non-parallel |
| 3 | ~1887-1898 | `total_2d` | XZ | Parallel (PE0 gather) |
| 4 | ~1940-1950 | `local_2d` | XZ | Non-parallel |
| 5 | ~2132-2142 | `total_2d` | YZ | Parallel (PE0 gather) |
| 6 | ~2185-2195 | `local_2d` | YZ | Non-parallel |
| 7 | ~2231-2245 | `local_2d_sections` | XY | Parallel netCDF4 |
| 8 | ~2288-2302 | `local_2d_sections` | XZ | Parallel netCDF4 |
| 9 | ~2352-2365 | `local_2d_sections` | YZ | Parallel netCDF4 |

**Pattern** (same for all 9):
```fortran
nan_count = COUNT( <array> /= <array> ) + COUNT( <array> > 1.0E30_wp ) + &
             COUNT( <array> < -1.0E30_wp )
IF ( nan_count > 0 )  THEN
   WRITE( message_string, * ) 'data_output_2d: NaN/Inf detected in variable ',&
                              TRIM( do2d(av,ivar) ), ''' (', nan_count, &
                              ' cells), replaced with 0.0 for 2D output'
   CALL message( 'data_output_2d', 'PAC1000', 0, 1, 0, 6, 0 )
   WHERE( <array> /= ... ) <array> = 0.0_wp
ENDIF
```

---

## File 4: `data_output_profiles.f90`

**Purpose**: Writes vertical profile output (`data_output_pr`)

### Change 4.1 — Import `message_string`

**Line**: ~30 (in `USE control_parameters`)
**Why**: Required — without this, compilation fails with `error #6404`

**Original**:                                      | **Modified**:
---                                                | ---
```fortrath                                           | ```fortran
USE control_parameters,                             | USE control_parameters,
    ONLY: average_count_pr,                         |     ONLY: average_count_pr,
           averaging_interval_pr,                   |            averaging_interval_pr,
           coupling_start_time,                     |            coupling_start_time,
           dopr_n,                                  |            dopr_n,
           dopr_time_count,                         |            dopr_time_count,
           time_since_reference_point                |            message_string,
```                                                  |            time_since_reference_point
                                                    | ```

### Change 4.2 — Add `nan_count` variable

**Line**: ~60

### Change 4.3 — NaN protection (Initial profiles, `hom`)

**Lines**: ~147-160 | **Array**: `hom(nzb:nzt+1,1,dopr_initial_index(i),sr)`

### Change 4.4 — NaN protection (Time-averaged, `hom_sum`)

**Lines**: ~200-218 | **Array**: `hom_sum(nzb:nzt+1,dopr_index(i),sr)`

---

## File 5: `data_output_tseries.f90`

**Purpose**: Writes time series output (`dt_dots`)

### Change 5.1 — Import `message_string`

**Line**: ~28 (in `USE control_parameters`)

**Original**:                                      | **Modified**:
---                                                | ---
```fortran                                         | ```fortran
USE control_parameters,                             | USE control_parameters,
    ONLY: dots_time_count,                          |     ONLY: dots_time_count,
           time_since_reference_point                |            message_string,
```                                                  |            time_since_reference_point
                                                    | ```

### Change 5.2 — Add `nan_count` variable

**Line**: ~60

### Change 5.3 — NaN protection (Time series, `ts_value`)

**Lines**: ~103-113 | **Array**: `ts_value(i,sr)` — scalar value (not array)

**Note**: Uses scalar `IF` instead of `COUNT()` because `ts_value` is a single value.

```fortran
nan_count = 0
IF ( ts_value(i,sr) /= ts_value(i,sr) )  nan_count = nan_count + 1
IF ( ts_value(i,sr) > 1.0E30_wp )  nan_count = nan_count + 1
IF ( ts_value(i,sr) < -1.0E30_wp )  nan_count = nan_count + 1
IF ( nan_count > 0 )  THEN
   WRITE( message_string, * ) 'data_output_tseries: NaN/Inf detected', &
                              ' (', nan_count, ' cells), replaced for time series'
   CALL message( 'data_output_tseries', 'PAC1000', 0, 1, 0, 6, 0 )
   ts_value(i,sr) = 0.0_wp
ENDIF
```

---

## File 6: `chemistry_model_mod.f90`

**Purpose**: Gas-phase chemistry ODE integration and data output

### Change 6.1 — Array bounds fix

**Line**: ~1392
**Bug**: `mask_size(mid,3)` → should be `mask_size_l(mid,3)`

**Original**:                           | **Modified**:
---                                     | ---
`DO k = 1, mask_size(mid,3)` ← BUG      | `DO k = 1, mask_size_l(mid,3)` ← FIXED

**Why**: Used global mask size instead of local per-processor size in `chem_data_output_mask`. In PALM's 2D decomposition (x,y only), z is not decomposed, so both are equal — but semantically incorrect.

### Change 6.2 — Rosenbrock floor (vector path)

**Lines**: ~2249-2254
**Subroutine**: `chem_integrate` (vectorized)

**Original**: No concentration floor

**Modified**:
```fortran
!-- Apply minimum concentration floor of 1.0 molecule/cm3 to ALL species
!-- to prevent Rosenbrock integrator failure. WRF-Chem input contains
!-- near-zero values (down to 1e-30 ppmv) for some species at certain
!-- grid points, making the ODE Jacobian near-singular.
!-- 1.0 molecule/cm3 is negligible for atmospheric chemistry.
                WHERE( tmp_conc(ks:ke,lsp) < 1.0_wp )  tmp_conc(ks:ke,lsp) = 1.0_wp
```

### Change 6.3 — Rosenbrock floor (scalar path)

**Lines**: ~2356-2361
**Subroutine**: `chem_integrate_ij` (point-wise)

Same fix as 6.2 for the scalar ODE path.

---

## File 7: `salsa_mod.f90`

**Purpose**: SALSA aerosol module

### Change 7.1 — NaN cleanup in dynamic driver

**Lines**: ~2288-2292
**Why**: Dynamic driver may contain NaN fill values in mass fractions

```fortran
          WHERE ( pmf2a /= pmf2a )  pmf2a = 0.0_wp
          WHERE ( pmf2b /= pmf2b )  pmf2b = 0.0_wp
          pmf2a = MAX( pmf2a, 0.0_wp )
          pmf2b = MAX( pmf2b, 0.0_wp )
!
!--       Store surface-level mass fractions for header output
          used_mass_fracs_a(:) = pmf2a(nzb+1,:)
          used_mass_fracs_b(:) = pmf2b(nzb+1,:)
```

### Change 7.2 — Add `used_mass_fracs_a/b` variables

**Lines**: ~473-479 | Module-level variables to store actual mass fractions in use

### Change 7.3 — Header: traffic conditional

**Lines**: ~1380-1384 | Only print traffic mass fractions when `is_traffic=.TRUE.`

```fortran
    IF ( is_traffic )  THEN
       WRITE( io, 151 ) mass_fracs_a_tra, mass_fracs_b_tra
    ENDIF
```

### Change 7.4 — Header: dynamic driver mass fractions

**Lines**: ~1374-1380 | Show actual values from dynamic driver

```fortran
    IF ( init_aerosol_type == 0 )  THEN
       WRITE( io, 15 ) ncc, listspec, mass_fracs_a, mass_fracs_b
    ELSE
       WRITE( io, 152 ) ncc, listspec, used_mass_fracs_a, used_mass_fracs_b
    ENDIF
```

### Change 7.5 — New format 152

**Lines**: ~1440-1444
```fortran
152 FORMAT (/'    Number of chemical components used: ', I1,/                             &
             '       Species: ',12(A6),/                                                   &
             '    Initial mass fractions (from dynamic driver, surface level):',/          &
             '       a-bins: ', 12(F10.7),/                                                &
             '       b-bins: ', 12(F10.7))
```

### Change 7.6 — Increased decimal places

**Format**: `F6.3` → `F10.7` in formats 15, 151, 152

---

### Change 7.7 — Non-convergence robustness in equilibration (subranges 1 & 2)

**Lines**: ~5088 (subrange 1a), ~5175 (subrange 2)
**Why**: During restart initialization, `salsa_driver` re-equilibrates the saved (evolved) aerosol state. At saturated RH (internally capped at 0.98 in `equilibration`), the ZSR/Kelvin fixed-point iteration can oscillate (e.g., bin 6). Previously the 1000-iteration guard called `message('SAL0029'/'SAL0030', 1, 2, 0, 6, 0)`: the message was only output on PE 0 (silent elsewhere), and `requested_action=1` triggered `local_stop` → collective `MPI_FINALIZE` on only the affected PEs while the rest kept running → MPI deadlock (restart "hang"). Now the loop breaks out with the last computed wet diameter and writes a `SALSA-WARN` line to the per-PE DEBUG file.

```fortran
             counti = counti + 1
             IF ( counti > 1000 )  THEN
                WRITE( 9, * )  'SALSA-WARN equilibration subrange 1 no convergence: bin=', ib, &
                               ' rh=', zrh, ' dwet=', zdwet, ' (using last value)'
                EXIT
             ENDIF
```

### Change 7.8 — Non-convergence robustness in distr_update

**Lines**: ~8987
**Why**: The particle bin-redistribution loop can fail to converge within 100 iterations on certain evolved states. Previously `SAL0037` deadlocked exactly as described in Change 7.7. Now it writes a `SALSA-WARN` line and exits, keeping the current bin state.

```fortran
       counti = counti + 1
       IF ( counti > 100 )  THEN
          WRITE( 9, * )  'SALSA-WARN distr_update not converged after 100 iterations (using current bin state)'
          EXIT
       ENDIF
```

### Change 7.9 — Clean visible abort for genuine distr_update error (SAL0036)

**Lines**: ~8984
**Why**: The genuine `zvfrac < 0` error guard previously used `requested_action=1` (silent subset-abort → deadlock). Changed to `requested_action=2` (MPI_ABORT on comm2d), `output_on_pe=-1` (message on all PEs), `file_id=9` (per-PE DEBUG file), `flush_file=1` — so if it fires, the whole job aborts cleanly with the message visible in every PE's DEBUG file.

**Original**: `CALL message( 'salsa_mod: distr_update', 'SAL0036', 1, 2, 0, 6, 0 )`
**Modified**: `CALL message( 'salsa_mod: distr_update', 'SAL0036', 2, 2, -1, 9, 1 )`

---

### Change 7.10 — Floating-invalid guard in `set_nest_mass` (offline nesting, Regime 2a)

**Lines**: ~18150-18190 (USER_CODE)
**Why**: Intermittent restart crashes with `forrtl: error (65): floating invalid` at `salsa_mod.f90:18046` (`set_nest_mass`) during offline-nesting init. The statement
`to_mconc = MAX(0.0_wp, mf2a/(SUM(mf2a)+eps)) * to_nconc * core * prho`
can raise a NaN via `0/0`, `Inf/Inf` (division) or `0*Inf` (multiplication). The `+ eps` only guards an exact `0.0` denominator, not NaN/Inf.

**Modified**: guard checks EVERY operand (`mf2a`, `sumf`, `to_nconc`, `core`) for NaN (`x /= x`) and ±Inf (`ABS(x) > HUGE(1.0_wp)`), plus `sumf <= 0`. On any bad cell: `to_mconc = 0.0_wp` + log line (species, boundary, grid `i/k/t`, bin, values; capped at 30 lines + per-species summary). The `ELSE` branch is identical to the original, so valid data is untouched.

```fortran
sumf = SUM( mf2a(t,k,i,:) )
IF ( .NOT. ( mf2a(t,k,i,ic)         == mf2a(t,k,i,ic) )         .OR.  &
     .NOT. ( sumf == sumf )                                    .OR.  &
     .NOT. ( to_nconc(t,k,i,ib)     == to_nconc(t,k,i,ib) )    .OR.  &
     ABS( mf2a(t,k,i,ic) )          > HUGE( 1.0_wp )           .OR.  &
     ABS( sumf )                    > HUGE( 1.0_wp )           .OR.  &
     ABS( to_nconc(t,k,i,ib) )      > HUGE( 1.0_wp )           .OR.  &
     ABS( core(ib) )                > HUGE( 1.0_wp )           .OR.  &
     sumf <= 0.0_wp )  THEN
   to_mconc(t,k,i,icc) = 0.0_wp
   ! ... log: species / boundary / i / k / t / ibin + mf2a / sumf / nconc / core / prho ...
ELSE
   to_mconc(t,k,i,icc) = MAX( 0.0_wp, mf2a(t,k,i,ic) / ( sumf + eps ) ) *   &
                         to_nconc(t,k,i,ib) * core(ib) * prho
ENDIF
```

### Change 7.11 — Floating-invalid guard in Regime 1 (`pmf1a` divisions)

**Lines**: ~17912-17994 (SO4 block ~17912-17951, OC block ~17954-17994)
**Why**: The same class of division (`pmf1a = mf2a / (mf2a[1] + mf2a[2] + eps)`) in `nesting_offl_aero_mass`, executed just before `set_nest_mass`, could equally crash. Added the identical guard on `mf2a(cc1)`, `mf2a(cc2)`, `den1 = mf2a(cc1)+mf2a(cc2)+eps`, `to_nconc`, `core`, and `den1 <= 0`.

**Modified**: bad cell → `to_mconc = 0.0_wp` + log (species SO4/OC, boundary, i/k/t, bin, values); `ELSE` = original arithmetic. Added declarations `den1`, `n_bad1`, `n_bad1_log`, `bname1` (~lines 17751-17754), boundary-name setup (~17892-17902), and per-call summary (~18030-18036).

```fortran
den1 = mf2a(t,k,i,cc_i2m(1)) + mf2a(t,k,i,cc_i2m(2)) + eps
IF ( .NOT. ( mf2a(t,k,i,cc_i2m(1)) == mf2a(t,k,i,cc_i2m(1)) ) .OR.  &
     .NOT. ( mf2a(t,k,i,cc_i2m(2)) == mf2a(t,k,i,cc_i2m(2)) ) .OR.  &
     .NOT. ( den1 == den1 )                                     .OR.  &
     .NOT. ( to_nconc(t,k,i,ib) == to_nconc(t,k,i,ib) )        .OR.  &
     ABS( mf2a(t,k,i,cc_i2m(1)) ) > HUGE( 1.0_wp )             .OR.  &
     ABS( mf2a(t,k,i,cc_i2m(2)) ) > HUGE( 1.0_wp )             .OR.  &
     ABS( den1 ) > HUGE( 1.0_wp )                              .OR.  &
     ABS( to_nconc(t,k,i,ib) ) > HUGE( 1.0_wp )                .OR.  &
     ABS( core(ib) ) > HUGE( 1.0_wp )                          .OR.  &
     den1 <= 0.0_wp )  THEN
   to_mconc(t,k,i,ic) = 0.0_wp
   ! ... log (SO4 or OC; boundary; i/k/t; ibin; values) ...
ELSE
   pmf1a = mf2a(t,k,i,cc_i2m(1)) / den1
   to_mconc(t,k,i,ic) = pmf1a * to_nconc(t,k,i,ib) * core(ib) * arhoh2so4
ENDIF
```

### Change 7.12 — Effect on simulation & output

- **Valid data**: no effect — the `ELSE` branch computes identical values (bit-identical simulation/output).
- **Invalid cell**: replaced by `0.0` (physical lower bound) + log; prevents NaN propagation; run continues instead of aborting.
- **Does not fix the root cause** — the log lines (species / boundary / time-slice / values) are the diagnostic to locate the bad value in the driving or restart data.
- **Activation**: requires rebuilding the job's executable.

---

## File 8: `restart_data_mpi_io_mod.f90`

**Purpose**: MPI-I/O for restart files (save/load of model state)

### Change 8.1 — Increase `max_nr_arrays`

**Line**: ~238
**Why**: Dual-population SALSA (traffic + non-traffic) doubles the number of aerosol restart arrays. The original limit of 600 was exceeded during restart writes, causing out-of-bounds access in `array_names`/`array_offset`. Increased to 5000.

**Original**: `INTEGER, PARAMETER :: max_nr_arrays = 600`
**Modified**: `INTEGER, PARAMETER :: max_nr_arrays = 5000`

---

## File 9: `pmc_interface_mod.f90`

**Purpose**: Domain nesting interface routines (parent-child nesting, interpolation/anterpolation)

### Change 9.1 — SALSA nesting array dimension (dual-population)

**Lines**: ~3047-3048 (`pmci_set_array_pointer`)
**Why**: Dual-population SALSA uses 5D nesting arrays with an explicit population index. Index 1 = non-traffic population; traffic variables are handled separately.

**Original**: `p_3d_sec => nconc_2(:,:,:,n)` / `p_3d_sec => mconc_2(:,:,:,n)`
**Modified**: `p_3d_sec => nconc_2(:,:,:,1,n)` / `p_3d_sec => mconc_2(:,:,:,1,n)`

---

## Files to Copy to HPC

```
JOBS/salsa_tra/USER_CODE/
├── data_output_mask.f90         — NaN protection for masked output
├── data_output_3d.f90           — NaN protection for 3D output
├── data_output_2d.f90           — NaN protection for 2D sections
├── data_output_profiles.f90     — NaN protection for profiles
├── data_output_tseries.f90      — NaN protection for time series
├── chemistry_model_mod.f90      — Rosenbrock floor + mask_size fix
├── salsa_mod.f90                — NaN cleanup + header improvements + restart non-convergence robustness
├── restart_data_mpi_io_mod.f90  — max_nr_arrays 600 → 5000
└── pmc_interface_mod.f90        — dual-population nesting array dimensions
```
