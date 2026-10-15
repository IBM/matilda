"""Background execution of harmonic persistence computations.

Each job runs :func:`matilda.api.service.compute_harmonic_persistence` in its own
subprocess, so a job that crashes, runs out of memory or exceeds its time limit fails
on its own without affecting the server or other jobs. Jobs and their results are kept
in memory.
"""

import multiprocessing
import resource
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

from .service import compute_harmonic_persistence

QUEUED = "queued"
RUNNING = "running"
SUCCEEDED = "succeeded"
FAILED = "failed"
CANCELLED = "cancelled"

FINISHED = (SUCCEEDED, FAILED, CANCELLED)


def _run_in_subprocess(connection, params, max_memory_mb):
    """compute one job and send ("ok", result) or ("error", message) through `connection`"""
    try:
        if max_memory_mb is not None:
            limit = int(max_memory_mb) * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        connection.send(("ok", compute_harmonic_persistence(**params)))
    except MemoryError:
        connection.send(("error", "out of memory (limit {} MB)".format(max_memory_mb)))
    except Exception as e:
        connection.send(("error", "{}: {}".format(type(e).__name__, e)))
    finally:
        connection.close()


def _timestamp(t):
    return None if t is None else datetime.fromtimestamp(t, timezone.utc).isoformat()


class Job:
    """state of one submitted computation"""

    def __init__(self, params):
        self.id = uuid.uuid4().hex
        self.params = params
        self.status = QUEUED
        self.created_at = time.time()
        self.started_at = None
        self.finished_at = None
        self.error = None
        self.result = None
        self.process = None

    def describe(self):
        """JSON-ready status, with a summary of the result once the job has succeeded"""
        summary = None
        if self.result is not None:
            summary = {
                "input": self.result["input"],
                "n_simplices": self.result["n_simplices"],
                "n_bars": len(self.result["barcode"]),
                "n_representatives": len(self.result["harmonic"]),
                "timings": self.result["timings"],
            }
        return {
            "id": self.id,
            "status": self.status,
            "created_at": _timestamp(self.created_at),
            "started_at": _timestamp(self.started_at),
            "finished_at": _timestamp(self.finished_at),
            "error": self.error,
            "summary": summary,
        }


class JobManager:
    """runs jobs in subprocesses, at most `max_workers` at a time

    Parameters
    ----------
    max_workers : int, optional
        number of jobs running concurrently, by default 1
    job_timeout : float, optional
        seconds after which a running job is killed, by default no limit
    max_memory_mb : int, optional
        address-space limit of each job's subprocess, by default no limit
    max_finished_jobs : int, optional
        number of finished jobs kept; the oldest are forgotten first, by default 100
    start_method : str, optional
        multiprocessing start method, by default "spawn", which is safe to use from a
        threaded server
    """

    def __init__(
        self,
        max_workers=1,
        job_timeout=None,
        max_memory_mb=None,
        max_finished_jobs=100,
        start_method="spawn",
    ):
        self.job_timeout = job_timeout
        self.max_memory_mb = max_memory_mb
        self.max_finished_jobs = max_finished_jobs
        self._context = multiprocessing.get_context(start_method)
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="matilda-job"
        )
        self._jobs = {}
        self._lock = threading.Lock()

    def submit(self, params):
        """queue a job computing `compute_harmonic_persistence(**params)` and return it"""
        job = Job(params)
        with self._lock:
            self._jobs[job.id] = job
            self._forget_old_jobs()
        self._executor.submit(self._run, job)
        return job

    def get(self, job_id):
        """the job with this id, or None"""
        with self._lock:
            return self._jobs.get(job_id)

    def list_jobs(self):
        """all known jobs, oldest first"""
        with self._lock:
            return list(self._jobs.values())

    def cancel(self, job_id):
        """stop a queued or running job; return False if the job is unknown"""
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return False
            if job.status in (QUEUED, RUNNING):
                job.status = CANCELLED
                job.finished_at = time.time()
                if job.process is not None:
                    job.process.kill()
            return True

    def delete(self, job_id):
        """cancel a job and forget it; return False if the job is unknown"""
        if not self.cancel(job_id):
            return False
        with self._lock:
            self._jobs.pop(job_id, None)
        return True

    def shutdown(self):
        """cancel all unfinished jobs and stop the worker threads"""
        for job in self.list_jobs():
            self.cancel(job.id)
        self._executor.shutdown(wait=False, cancel_futures=True)

    def _forget_old_jobs(self):
        finished = [job for job in self._jobs.values() if job.status in FINISHED]
        for job in finished[: max(0, len(finished) - self.max_finished_jobs)]:
            del self._jobs[job.id]

    def _run(self, job):
        with self._lock:
            if job.status != QUEUED:
                return
            job.status = RUNNING
            job.started_at = time.time()

        receiver, sender = self._context.Pipe(duplex=False)
        process = self._context.Process(
            target=_run_in_subprocess,
            args=(sender, job.params, self.max_memory_mb),
            daemon=True,
        )
        process.start()
        sender.close()
        with self._lock:
            job.process = process
            if job.status == CANCELLED:
                process.kill()

        outcome = None
        try:
            if receiver.poll(self.job_timeout):
                outcome = receiver.recv()
            else:
                outcome = (
                    "error",
                    "the job exceeded the time limit of {} s".format(self.job_timeout),
                )
        except EOFError:
            pass
        finally:
            if process.is_alive():
                process.kill()
            process.join()
            receiver.close()

        if outcome is None:
            outcome = (
                "error",
                "the job's process exited unexpectedly (exit code {}), "
                "possibly out of memory".format(process.exitcode),
            )

        with self._lock:
            job.process = None
            if job.status == CANCELLED:
                return
            job.finished_at = time.time()
            if outcome[0] == "ok":
                job.status = SUCCEEDED
                job.result = outcome[1]
            else:
                job.status = FAILED
                job.error = outcome[1]
            self._forget_old_jobs()
