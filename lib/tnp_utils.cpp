//File with C++ utilities for T&P analysis

//function to avoid memory leakage in PyROOT from RooFit
void free_memory_RooFitResult(RooFitResult* ptr) {
  delete ptr;
}
