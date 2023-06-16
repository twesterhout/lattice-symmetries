module Utils {

private use List;
private use FFI;
import IO;

proc initRuntime() {
  coforall loc in Locales do on loc {
    ls_hs_init();
  }
}

proc deinitRuntime() {
  coforall loc in Locales do on loc {
    ls_hs_exit();
  }
}

class RoseTree {
  type eltType;
  var current : eltType;
  var children : list(shared RoseTree(eltType));

  proc init(x : ?eltType) {
    this.eltType = eltType;
    this.current = x;
  }
  proc init(x : ?eltType, ys) {
    this.eltType = eltType;
    this.current = x;
    this.children = ys;
  }

  proc addChild(x : shared RoseTree(eltType)) {
    this.children.pushBack(x);
  }
  proc addChild(x : eltType) {
    addChild(new shared RoseTree(x));
  }

  proc getChild(x : eltType) ref {
    for c in children {
      if c.current == x then return c;
    }
    halt("child " + x:string + " not found");
  }

  // proc map(fn : proc(_ : eltType) : ?t) : shared RoseTree(t) {
  //   var outTree = new shared RoseTree(fn(current));
  //   for c in children do
  //     outTree.addChild(c.map(fn));
  //   return outTree;
  // }

  proc writeChildren(f, indent : string) throws {
    if !children.isEmpty() {
      for (x, i) in zip(children, 0 ..) {
        const currentPrefix = if i == children.size - 1 then "└─ " else "├─ ";
        const childPrefix = if i == children.size - 1 then "   " else "│  ";
        f.write(indent, currentPrefix, x.current, '\n');
        x.writeChildren(f, indent + childPrefix);
      }
    }
  }
  override proc writeThis(f) throws {
    f.write(current, '\n');
    writeChildren(f, indent = "  ");
  }
}

record TimingResult {
  type eltType;
  var func : string;
  var stat : eltType;

  proc init(func : string, stat : ?eltType) {
    this.eltType = eltType;
    this.func = func;
    this.stat = stat;
  }

  proc writeThis(f) throws {
    f.write(func, ": ", stat);
  }
}

record TimingResultAccumulator {
  type eltType;
  var func : string;
  var _count : int;
  var _mean : eltType;
  var _M2 : eltType;

  proc init(func : string, x : ?eltType) {
    this.eltType = eltType;
    this.func = func;
    this._count = 1;
    this._mean = x;
    this._M2 = 0;
  }
  proc init(x : TimingResult(?eltType)) {
    init(x.func, x.stat);
  }

  proc update(x : TimingResult(eltType)) {
    if func != x.func then
      halt("cannot combine TimingResults with different func attributes: " + func + " != " + x.func);
    _count += 1;
    const delta = x.stat - _mean;
    _mean += delta / _count;
    const delta2 = x.stat - _mean;
    _M2 += delta * delta2;
  }

  inline proc mean { return if _count == 0 then (0.0 / 0.0):eltType else _mean; }
  inline proc std { return if _count < 2 then (0.0 / 0.0):eltType else _M2 / _count; }
}

type TimingTree = shared RoseTree(TimingResult(real));

proc timingTree(func, result, children) {
  return new shared RoseTree(
    new TimingResult(func, result),
    [(f, r) in children] new shared RoseTree(new TimingResult(f, r))
  );
}
proc timingTree(func, result) {
  return new shared RoseTree(new TimingResult(func, result));
}

proc meanAndErr(xs : [] real) {
  const mean = (+ reduce xs) / xs.size:real;
  const variance =
    (1.0 / xs.size:real)
      * (+ reduce ([i in xs.domain] (xs[i] - mean) * (xs[i] - mean)));
  const err = round(100 * sqrt(variance)) / 100;
  return (mean, err);
}

proc combineTimingTrees(ref dest : shared RoseTree(TimingResultAccumulator(?eltType)),
                        const ref tree : shared RoseTree(TimingResult(eltType))) {
  dest.current.update(tree.current);
  for (destChild, child) in zip(dest.children, tree.children) {
    combineTimingTrees(destChild, child);
  }
}

proc timingTreeToAccumulator(tree : TimingTree) : shared RoseTree(TimingResultAccumulator(real)) {
  var outTree = new shared RoseTree(new TimingResultAccumulator(tree.current));
  for c in tree.children do
    outTree.addChild(timingTreeToAccumulator(c));
  return outTree;
}

proc timingTreeFromAccumulator(tree : shared RoseTree(TimingResultAccumulator(real))) : TimingTree {
  var outTree = new shared RoseTree(new TimingResult(tree.current.func, tree.current.mean));
  for c in tree.children do
    outTree.addChild(timingTreeFromAccumulator(c));
  return outTree;
}

proc combineTimingTrees(trees : [] shared RoseTree(TimingResult(real))) {
  type eltType = real;
  var arr = trees;
  var dest = timingTreeToAccumulator(arr[0]);
  // var dest = arr[0].map(
  //   proc(x : TimingResult(eltType)) : TimingResultAccumulator(eltType) {
  //     return new TimingResultAccumulator(x);
  //   });
  for tree in arr[1..] do
    combineTimingTrees(dest, tree);
  return timingTreeFromAccumulator(dest);
  // return dest.map(
  //   proc(x : TimingResultAccumulator(eltType)) : TimingResult(eltType) {
  //     return new TimingResult(x.func, x.mean);
  //   });
}


} // module Utils
