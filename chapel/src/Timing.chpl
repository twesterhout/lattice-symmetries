module Timing {

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
    assert(r.count > 0);
    this.totalTime += r.time;
    this.minTime = min(this.minTime, r.time);
    this.maxTime = max(this.maxTime, r.time);
    this.totalCount += r.count;
    this.minCount = min(this.minCount, r.count);
    this.maxCount = max(this.maxCount, r.count);
    this._n += 1;
  }
}

proc TimingCombined.serialize(writer : fileWriter(?), ref serializer) throws {
  var ser = serializer.startRecord(writer, "TimingCombined", 9);
  try ser.writeField("totalTime", totalTime);
  try ser.writeField("avgTime", totalTime / _n);
  try ser.writeField("minTime", minTime);
  try ser.writeField("maxTime", maxTime);
  try ser.writeField("totalCount", totalCount);
  try ser.writeField("avgCount", totalCount:real / _n);
  try ser.writeField("minCount", minCount);
  try ser.writeField("maxCount", maxCount);
  try ser.writeField("numTasks", _n);
  ser.endClass();
}

record TaskLocalStore {
  var store : map(string, TimingRecord, parSafe=false);

  proc ref add(param func, time : real, count : int = 1) {
    store._enter(); defer store._leave();

    var (found, slot) = store.table.findAvailableSlot(func);
    if found {
      ref value = store.table.table[slot].val;
      value += new TimingRecord(time, count);
    }
    else {
      store.table.fillSlot(slot, func, new TimingRecord(time, count));
    }
  }
}

class LocaleLocalStore {
  var store : map(chpl_taskID_t, TaskLocalStore, parSafe=true);

  proc add(param func, time : real) {
    const taskId = chpl_task_getId();

    store._enter(); defer store._leave();
    var (found, slot) = store.table.findAvailableSlot(taskId);
    if !found {
      store.table.fillSlot(slot, taskId, new TaskLocalStore());
    }
    ref value = store.table.table[slot].val;
    value.add(func, time);
  }
}

proc initLocaleLocalStore() {
  var obj = new unmanaged LocaleLocalStore();
  c_timing_set_locale_data(c_ptrTo(obj));
}

proc getLocaleLocalStore() {
  var obj : unmanaged LocaleLocalStore = (c_timing_get_locale_data():unmanaged LocaleLocalStore?)!;
  return obj;
}

record Recorder {
  param func;
  var timer : stopwatch;

  proc init(param func) {
    this.func = func;
    init this;
    timer.start();
  }

  proc ref deinit() {
    timer.stop();
    getLocaleLocalStore().add(func, timer.elapsed());
  }
}

proc recordTime(param func) { return new Recorder(func); };

proc summarize() {
  var r : map(string, TimingCombined);

  for table in getLocaleLocalStore().store.values() {
    for funcName in table.store.keys() {
      const ref timing = try! table.store[funcName];
      if r.contains(funcName) {
        r[funcName].combine(timing);
      }
      else {
        var t : TimingCombined;
        t.combine(timing);
        r.add(funcName, t);
      }
    }
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
