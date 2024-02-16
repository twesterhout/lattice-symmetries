module Timing {

private import Reflection;
private use CTypes;
private use IO;
private use Map;
private use Time;
private use Sort;

require "Timing.h";

extern proc chpl_task_getId(): chpl_taskID_t;
extern proc c_timing_get_locale_data(): c_ptr(void);
extern proc c_timing_set_locale_data(ptr : c_ptr(void));

record TimingRecord {
  var time : real;
  var count : int;
}

operator TimingRecord.+=(ref x : TimingRecord, t : real) {
  x.time += t;
  x.count += 1;
}
operator TimingRecord.+=(ref x : TimingRecord, const ref y : TimingRecord) {
  x.time += y.time;
  x.count += y.count;
}

record TimingCombined : writeSerializable {
  var totalTime : real;
  var minTime : real;
  var maxTime : real;

  var totalCount : int;
  var minCount : int;
  var maxCount : int;

  var _n : int;

  proc init() {
    this.totalTime = 0;
    this.minTime = inf;
    this.maxTime = 0;
    this.totalCount = 0;
    this.minCount = max(int);
    this.maxCount = 0;
    this._n = 0;
  }

  proc ref combine(const ref r : TimingRecord) {
    if r.count > 0 {
      this.totalTime += r.time;
      this.minTime = min(this.minTime, r.time);
      this.maxTime = max(this.maxTime, r.time);
      this.totalCount += r.count;
      this.minCount = min(this.minCount, r.count);
      this.maxCount = max(this.maxCount, r.count);
      this._n += 1;
    }
  }
}

proc TimingCombined.serialize(writer : fileWriter(?), ref serializer) throws {
  var ser = serializer.startRecord(writer, "TimingCombined", 9);
  try ser.writeField("totalTime", totalTime);
  try ser.writeField("avgTime", if _n > 0 then totalTime / _n else -1);
  try ser.writeField("minTime", if _n > 0 then minTime else -1);
  try ser.writeField("maxTime", if _n > 0 then maxTime else -1);
  try ser.writeField("totalCount", totalCount);
  try ser.writeField("avgCount", if _n > 0 then totalCount:real / _n else -1);
  try ser.writeField("minCount", if _n > 0 then minCount else -1);
  try ser.writeField("maxCount", if _n > 0 then maxCount else -1);
  try ser.writeField("numTasks", _n);
  ser.endClass();
}

// This part is pretty ugly, but at least the runtime overhead of accessing a
// field is much lower than doing a hash table lookup 😄
record Store {
  var _enumerateStatesProjected : TimingRecord;
  var _enumerateStatesUnprojected : TimingRecord;
  var basisStatesToIndices : TimingRecord;
  var computeOffDiagNoProjection : TimingRecord;
  var computeOffDiagWithProjection : TimingRecord;
  var convertOffDiagToCsr : TimingRecord;
  var determineEnumerationRanges : TimingRecord;
  var extractDiag : TimingRecord;
  var isRepresentative : TimingRecord;
  var localProcessExperimental : TimingRecord;
  var perLocaleDiagonal : TimingRecord;
  var perLocaleMatrixVector : TimingRecord;
  var perLocaleOffDiagonal : TimingRecord;
}


record TaskLocalStore {
  var store : Store;

  inline proc ref add(param func, time : real, count : int = 1) {
    ref field = Reflection.getFieldRef(this.store, func);
    field += new TimingRecord(time, count);
  }
}

class LocaleLocalStore {
  var dom : domain(1);
  var arr : [dom] TaskLocalStore;

  proc init(capacity : int = 10240) {
    this.dom = {0 ..# capacity};
  }

  inline proc add(param func, time : real) {
    const taskId = chpl_task_getId():int;
    if taskId >= dom.size then
      halt(try! "taskId: {}, but capacity={}".format(taskId, dom.size));
    arr[taskId].add(func, time);
  }
}

proc initLocaleLocalStore() {
  var obj = new unmanaged LocaleLocalStore();
  c_timing_set_locale_data(c_ptrTo(obj));
}

inline proc getLocaleLocalStore() {
  var obj : unmanaged LocaleLocalStore = (c_timing_get_locale_data():unmanaged LocaleLocalStore?)!;
  return obj;
}

record Recorder {
  param func;
  var timer : stopwatch;

  inline proc init(param func) {
    this.func = func;
    init this;
    timer.start();
  }

  inline proc ref deinit() {
    timer.stop();
    getLocaleLocalStore().add(func, timer.elapsed());
  }
}

config param kRecordTiming = true;

inline proc recordTime(param func) where kRecordTiming { return new Recorder(func); }
inline proc recordTime(param func) where !kRecordTiming { return nil; }

proc summarizeRecursive(ref r, const ref table, param k, param numFields) : none {
  if k >= numFields {
    return;
  }
  else {
    param funcName = Reflection.getFieldName(table.type, k);
    const ref timing = Reflection.getField(table, k);
    if r.contains(funcName) {
      r[funcName].combine(timing);
    }
    else {
      var t : TimingCombined;
      t.combine(timing);
      r.add(funcName, t);
    }
    summarizeRecursive(r, table, k + 1, numFields);
  }
}


proc summarize() {
  var r : map(string, TimingCombined);

  const ref localeStore = getLocaleLocalStore();
  for x in localeStore.arr {
    const ref table = x.store;
    param numFields = Reflection.getNumFields(table.type);
    summarizeRecursive(r, table, 0, numFields);
  }

  var arr = r.toArray();

  record Comparator { proc key(x) : real { return -x[1].totalTime; } }
  sort(arr, comparator=new Comparator());
  return arr;
}


// Initialize on every locale
coforall loc in Locales do on loc do
  initLocaleLocalStore();

}

/*
module Main {
  use Timing;
  use IO, JSON;
  use Reflection;

  proc foo() {
    var _t = recordTime(getRoutineName());
    // if (chpl_taskID_t:int % 2 == 0) {
    writeln("foo");
    // }
  }

  proc main() {
    writeln("hello world!");

    forall i in 1..100 {
      foo();
    }

    writeln(getLocaleLocalStore());
    try! stdout.withSerializer(serializer=new jsonSerializer()).writeln(summarize());
  }
}
*/
