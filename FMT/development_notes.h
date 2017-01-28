/* Notes about FMT to help me keep track of what's going on.

 ******************
 Operational Status
 ******************

 Top priority item is to move to a kinodynamic version.
 
 Planning seems _ok_, but there's a big issue where the paths never seem to be executed all the way.
 
 Hitting limit cycles within the planned path. Move base issue?
 
 
 
 **************
 Current focus:
 **************
 
 Make explorer recompute only if needed (e.g. if LOS extension of existing path gets you to the goal, then skip re-planning)
 
 FMT sometimes remembers the previous goal (side effect of not resetting neighborhoods correctly). I think I've caught these now.
 
 
 
 

 *****************
 Summary of Issues
 *****************
 
 Polynomial issues:
 > Cell checking is wasteful
 > Cost of the polynomial is not computed correctly
 > Should write a suite of basic tests to see if things function the way I expect them to
 
 FMT issues:
 > Collision checking should be split into cached (for pre-generated tree) and on-line (for start/end connections). Cacheing the on-line portion is actually more expensive.
 > Need a rewiring step to smooth the output of the tree
 
 Explorer issues:
 > Seeing switching behavior again
 > Consider dynamic allocation for items with a large footprint (e.g. polyPath)
 
 Simulator issues:
 > Move_base is a big hack.
 > Consider a propogate dynamics function, so that we can decouple planning frequency from physics fidelity
 
 ********
 Comments
 ********
 
 (POLY) Timing is _extremely_ sensitive to the dt parameter in the path checking. See other note about not cacheing the on-line paths.
 
 (POLY) Look into `garbage value' claims by the analyzer. I think everything is ok, but check anyway.

 (POLY) Check that -1 is consistenly used as a INVALID_POLYNOMIAL cost, rather than MAXFLOAT

 (POLY) cheat in the distance returned by fit_polynomial right now. Should be done better. Same goes for the indices.
 
 (POLY) Reverse flag doesn't make sense given the initial conditions, should be eliminated

 (POLY/UTILS) Should make polystate a class with a copy constructor? Right now just have a copy function.
 
 
 (UTILS) Use appropriate classes to replace "utilities.h", e.g. map, graph, search, etc.
 Smoother interface - should make polynomial agnostic so we can easily switch to a 2pbvp approach/exact approach.
 
 (FMT) Downside of current FMT_Points structure - must represent x,y,z coordinates in FMT_ORD derivatives. Couldhave a more compact representation where we limit the order unevenly between x,y,z

 (FMT) Important to keep neighborhood lists sorted by index (rather than distance) to get speed-up in filtering operations (which occur very frequently).
 

 (FMT/POLY) Reconsider global MAX_POLY_CELLS in favor of dynamic allocation. Current practice leads to a lot of memory footprint.
 
 (POLY/EXP) Should redefine polypath vectors to be waaaay smaller.
 
 (SIM/EXP) Delta time is defined in both the explorer and simulator.
 
 (ALL) Reduce memory footprint:
 > MAX_POLY_CELLS is too small to be conservative (there are paths which use all of it), but at the same time the average path only uses half of the bytes. This wastes 65k for the small graph and 1.6M for the larger one. Significant, but not huge
 > Each PolyState object takes 520 bytes, and the current_plan container only uses ~2-5 of 1024 allocated objects. This is a 500k waste.
 > F_dense is around 85 MB, F_sparse is around 35MB accounting for the 2/3 of the memory allocated.
 
*/
