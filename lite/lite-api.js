// DrillPyze Lite: routes the front end's API calls to the browser engine in a Web Worker
// instead of the Python server.  Loaded before app.js, which checks for window.DRILLPYZE_LITE.
(function () {
  const worker = new Worker('engine/worker.js', { type: 'module' });
  const pending = new Map();
  let nextId = 1;
  worker.onmessage = (e) => {
    const { id, ok, data, error } = e.data, p = pending.get(id);
    if (!p) return;
    pending.delete(id);
    if (ok) p.resolve(data); else p.reject(new Error(error));
  };
  worker.onerror = (e) => {
    for (const p of pending.values()) p.reject(new Error(e.message || 'The simulation engine failed to start.'));
    pending.clear();
  };
  function call(path, body) {
    return new Promise((resolve, reject) => {
      const id = nextId++;
      pending.set(id, { resolve, reject });
      worker.postMessage({ id, path, body });
    });
  }
  window.DRILLPYZE_LITE = {
    api: call,
    // results download: the CSV is built in the worker and saved from a blob
    async download(session) {
      const { csv } = await call(`/api/results.csv?session=${session}`);
      const url = URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
      const a = Object.assign(document.createElement('a'), { href: url, download: 'drillpyze_results.csv' });
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    },
  };
})();
