/* Notes about FMT to help me keep track of what's going on.

 ******************
 Operational Status
 ******************
 Code runs. Paths returned are crap. Collision checking appears to fail most of the time.
 
 *****************
 Summary of Issues
 *****************
 
 Polynomial issues:
 > Wasteful checking of cells (get trevor's code and see if that is better)
 > Cost of the polynomial is not computed correctly
 > Should write a suite of basic tests to see if things function the way I expect them to
 
 FMT issues:
 > Collision checking does not seem to do anything
 > Polynomial paths returned as `good' do not seem good at all. Might be adding an extra segment to the visualization
 > Planner does not compensate for drift
 
 Explorer issues:
 > Seeing switching behavior again
 > Huge waste of memory in current implementation
 > Consider dynamic allocation for items with a large footprint (e.g. polyPath)
 
 Simulator issues:
 > Replanning occurs faster than the dynamics, so continuity constraints turn into inaction constraints.
 
 **************
 Current focus:
 **************
 
    Need to verify that the right cells are checked. Seems like collision checking is largely irrelevant to the paths chosen.
 
    Really should look into the drift version of FMT
 
 ********
 Comments
 ********
 (ALL) Reduce memory footprint
 
 (UTILS) Use appropriate classes to replace "utilities.h", e.g. map, graph, search, etc.
 Smoother interface - should make polynomial agnostic so we can easily switch to a 2pbvp approach/exact approach.
 
 (FMT) Downside of current FMT_Points structure - must represent x,y,z coordinates in FMT_ORD derivatives. Couldhave a more compact representation where we limit the order unevenly between x,y,z

 (FMT) Important to keep neighborhood lists sorted by index (rather than distance) to get speed-up in filtering operations (which occur very frequently).
 
 (POLY) cheat in the distance returned by fit_polynomial right now. Should be done better. Same goes for the indices.
 
 (POLY/UTILS) Should make polystate a class with a copy constructor? Right now just have a copy function.

 (FMT/POLY) Reconsider global MAX_POLY_CELLS in favor of dynamic allocation. Current practice leads to a lot of memory footprint.
 
 (FMT) Check whether there is a memory leak in reset_neighborhood. Might be an issue with the cells list by `delete by forgetting' approach.
 
 (POLY) Still something funky with path costs, but it seems like polyfmt is sort of working now.
 
 (POLY/EXP) Should redefine polypath vectors to be waaaay smaller.
 
 (SIM/EXP) Delta time is defined in both the explorer and simulator.
 
 (POLY) Reverse flag doesn't make sense given the initial conditions

*/