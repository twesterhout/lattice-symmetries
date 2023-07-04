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

record VarianceAccumulator {
  type eltType;
  var _count : int;
  var _mean : eltType;
  var _M2 : eltType;

  proc init(x : ?eltType) {
    this.eltType = eltType;
    this._count = 1;
    this._mean = x;
    this._M2 = 0;
  }

  proc update(x : eltType) {
    _count += 1;
    const delta = x - _mean;
    _mean += delta / _count;
    const delta2 = x - _mean;
    _M2 += delta * delta2;
  }

  inline proc mean { return if _count == 0 then (0.0 / 0.0):eltType else _mean; }
  inline proc std { return if _count < 2 then (0.0 / 0.0):eltType else _M2 / _count; }

  proc writeThis(f) throws {
    f.write(mean, " ± ", std);
  }
}


class RoseTree {
  type eltType;
  var func : string;
  var stat : eltType;
  var children : list(shared RoseTree(eltType));

  proc init(func : string, stat : ?eltType) {
    this.eltType = eltType;
    this.func = func;
    this.stat = stat;
  }
  proc init(func : string, stat : ?eltType, children) {
    this.eltType = eltType;
    this.func = func;
    this.stat = stat;
    this.children = children;
  }

  proc addChild(x : shared RoseTree(eltType)) {
    this.children.pushBack(x);
  }
  proc addChild(func : string, stat : eltType) {
    addChild(new shared RoseTree(func, stat));
  }

  proc getChild(func : string) ref {
    for c in children {
      if c.func == func then return c;
    }
    halt("child " + func:string + " not found");
  }

  proc toAccumulatorTree() : shared RoseTree(VarianceAccumulator(eltType)) {
    var tree = new shared RoseTree(func, new VarianceAccumulator(stat));
    for c in children do
      tree.addChild(c.toAccumulatorTree());
    return tree;
  }

  proc writeChildren(f, indent : string) throws {
    if !children.isEmpty() {
      for (x, i) in zip(children, 0 ..) {
        const currentPrefix = if i == children.size - 1 then "└─ " else "├─ ";
        const childPrefix = if i == children.size - 1 then "   " else "│  ";
        f.write(indent, currentPrefix, x.func, ": ", x.stat, '\n');
        x.writeChildren(f, indent + childPrefix);
      }
    }
  }
  override proc writeThis(f) throws {
    f.write(func, ": ", stat, '\n');
    writeChildren(f, indent = "  ");
  }
}


type TimingTree = shared RoseTree(real);

proc timingTree(func, result, children) {
  return new shared RoseTree(func, result,
    [(f, r) in children] new shared RoseTree(f, r));
}
proc timingTree(func, result) {
  return new shared RoseTree(func, result);
}

proc meanAndErr(xs : [] real) {
  const mean = (+ reduce xs) / xs.size:real;
  const variance =
    (1.0 / xs.size:real)
      * (+ reduce ([i in xs.domain] (xs[i] - mean) * (xs[i] - mean)));
  const err = round(100 * sqrt(variance)) / 100;
  return (mean, err);
}

proc fromAccumulatorTree(tree) : TimingTree {
  var outTree = new shared RoseTree(tree.func, tree.stat.mean);
  for c in tree.children do
    outTree.addChild(fromAccumulatorTree(c));
  return outTree;
}

proc combineTimingTrees(ref dest, const ref tree) {
  dest.stat.update(tree.stat);
  for (destChild, child) in zip(dest.children, tree.children) {
    combineTimingTrees(destChild, child);
  }
}
proc combineTimingTrees(arr : [] TimingTree) {
  if arr.isEmpty() then
    halt("trying to combine an empty list of timings");
  var dest = arr[0].toAccumulatorTree();
  for tree in arr[1..] do
    combineTimingTrees(dest, tree);
  return fromAccumulatorTree(dest);
}

} // module Utils
