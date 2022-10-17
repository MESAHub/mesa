#! Include : hdf5io
#! Purpose : support macros for hdf5io module
#!
#! Copyright 2021 Rich Townsend
#!
#! This file is part of the ForUM (Fortran Utility Modules)
#! package. ForUM is free software: you can redistribute it and/or
#! modify it under the terms of the GNU General Public License as
#! published by the Free Software Foundation, version 3.
#!
#! ForUM is distributed in the hope that it will be useful, but WITHOUT
#! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
#! License for more details.
#!
#! You should have received a copy of the GNU General Public License
#! along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!!
#! Standard macro to call HDF5
#!!

#:def HDF5_CALL(proc, *vars_pos, **vars_kw)
#:set vars = vars_pos + ['hdf_err'] + ['{:s}={:s}'.format(key, value) for key, value in vars_kw.items()]
#:set arg_list = '(' + ','.join(vars) + ')'
call ${proc}$${arg_list}$
if (hdf_err == -1) then
   call h5eprint_f(hdf_err)
   @:ABORT(error in call to ${proc}$)
endif
#:enddef
