! star_data.inc

include "star_data_step_input.inc"
! input for taking a step.
! modified during the step so that at end holds the new model.
! can be saved in models and photos.

include "star_data_step_work.inc"
! working storage for data derived from inputs while taking a step.
! can be undefined when start to take a step (e.g. by calling fill_arrays_with_NaNs).
! not saved in models or in photos.

include "star_data_procedures.inc"
! procedure pointers for hooks and such
