/***************************************************************************
                               morse_lj_dsf.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for acceleration of the lj/cut/coul/dsf pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 7/12/2012
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_MORSE_LJ_DSF_H
#define LAL_MORSE_LJ_DSF_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class MORLJDSF : public BaseCharge<numtyp, acctyp> {
 public:
  MORLJDSF();
  ~MORLJDSF();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    *
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, double **host_cutsq, double **host_lj1,
           double **host_lj2, double **host_lj3, double **host_lj4,
           double **host_offset, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen, double **host_cut_ljsq,
           const double host_cut_coulsq, double *host_special_coul,
           const double qqrd2e, const double e_shift, const double f_shift,
           const double alpha,double **host_morse1, double **host_r0,
           double **host_beta, double **host_d0);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = cutsq_vdw, lj1.w = cutsq
  UCL_D_Vec<numtyp4> lj1;
  /// lj3.x = lj3, lj3.y = lj4, lj3.z = offset
  UCL_D_Vec<numtyp4> lj3;
  /// mor1.x = cutsq, mor1.y = morse1, mor1.z = r0, mor1.w = beta
  UCL_D_Vec<numtyp4> mor1;
  /// mor2.x = d0, mor2.y = offset
  UCL_D_Vec<numtyp2> mor2;
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e;

 private:
  bool _allocated;
  numtyp _e_shift, _f_shift, _alpha, _cut_coulsq;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
