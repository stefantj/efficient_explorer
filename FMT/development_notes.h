/* Notes about FMT to help me keep track of what's going on.

 ******************
 Operational Status
 ******************
 Code runs. Paths returned are crap. Collision checking appears to fail most of the time. 
 
 Top priority item is to move to a drift sensitive version of FMT or an actual kinodynamic version.
 
 *****************
 Summary of Issues
 *****************
 
 Polynomial issues:
 > Wasteful checking of cells (get trevor's code and see if that is better)
 > Cost of the polynomial is not computed correctly
 > Should write a suite of basic tests to see if things function the way I expect them to
 
 FMT issues:
 > Collision checking does not seem to do anything
 > Might be adding an extra segment to the visualization
 > Planner does not compensate for drift
 > Collision checking should be split into cached (for pre-generated tree) and on-line (for start/end connections). Cacheing the on-line portion is actually more expensive.
 
 Explorer issues:
 > Seeing switching behavior again
 > Consider dynamic allocation for items with a large footprint (e.g. polyPath)
 
 Simulator issues:
 > Replanning occurs faster than the dynamics, so continuity constraints turn into inaction constraints.
 
 **************
 Current focus:
 **************
 
    Seems like collision checking is largely irrelevant to the paths chosen.
 
    Coefficient overflows are strange - look into this. Consider replacing the current stepping code with calls to `get_poly_der'?
 
    Really should look into the drift version of FMT


 
 ********
 Comments
 ********
 
 (POLY) Timing is _extremely_ sensitive to the dt parameter in the path checking. See other note about not cacheing the on-line paths.
 
 (POLY) Look into `garbage value' claims by the analyzer. I think everything is ok, but check anyway.
 
 (ALL) Reduce memory footprint:
 > MAX_POLY_CELLS is too small to be conservative (there are paths which use all of it), but at the same time the average path only uses half of the bytes. This wastes 65k for the small graph and 1.6M for the larger one. Significant, but not huge
 > Each PolyState object takes 520 bytes, and the current_plan container only uses ~2-5 of 1024 allocated objects. This is a 500k waste.
 > F_dense is around 85 MB, F_sparse is around 35MB accounting for the 2/3 of the memory allocated.

 
 (UTILS) Use appropriate classes to replace "utilities.h", e.g. map, graph, search, etc.
 Smoother interface - should make polynomial agnostic so we can easily switch to a 2pbvp approach/exact approach.
 
 (FMT) Downside of current FMT_Points structure - must represent x,y,z coordinates in FMT_ORD derivatives. Couldhave a more compact representation where we limit the order unevenly between x,y,z

 (FMT) Important to keep neighborhood lists sorted by index (rather than distance) to get speed-up in filtering operations (which occur very frequently).
 
 (POLY) cheat in the distance returned by fit_polynomial right now. Should be done better. Same goes for the indices.
 
 (POLY/UTILS) Should make polystate a class with a copy constructor? Right now just have a copy function.

 (FMT/POLY) Reconsider global MAX_POLY_CELLS in favor of dynamic allocation. Current practice leads to a lot of memory footprint.
 
 (POLY) Still something funky with path costs, but it seems like polyfmt is sort of working now.
 
 (POLY/EXP) Should redefine polypath vectors to be waaaay smaller.
 
 (SIM/EXP) Delta time is defined in both the explorer and simulator.
 
 (POLY) Reverse flag doesn't make sense given the initial conditions

*/
