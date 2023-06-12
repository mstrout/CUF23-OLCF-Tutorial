// stream-ep.chpl

// Example that can run stream on a cpu multicore and on any of the GPUs
// that node may have.
config var n = 1_000_000,
           alpha = 0.01;

coforall loc in Locales {
  on loc {
    cobegin {
      coforall gpu in here.gpus do on gpu {
        runTriad();
      }
      runTriad();        
    }
  }  
}
proc runTriad() {
   var A, B, C: [1..n] real;
   A = B + alpha * C;
}

